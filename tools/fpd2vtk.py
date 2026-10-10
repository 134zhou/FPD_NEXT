#!/usr/bin/env python3
"""
.fpd -> VTK 转换器。**离线运行，不占模拟的热路径。**

必须跑在 pvpython 下（需要 numpy + vtk）：
  PV=/home/doll/Software/ParaView-6.2.0-RC1-MPI-Linux-Python3.12-x86_64
  $PV/bin/pvpython tools/fpd2vtk.py out/run_*.fpd --config config/production.cfg -o vis/

产出：
  vis/<name>_<step>.vti   流体（ImageData，CellData: velocity 插值到格心）
  vis/<name>_<step>.vtp   粒子（PolyData，PointData: V / F / Ru）
  vis/<name>_<step>.vtm   多块，把上面两个绑在一起
  vis/<name>.pvd          整个时间序列 —— 【在 ParaView 里打开这一个文件即可】

可重复运行（边下载边转换）：产物齐全的文件直接跳过，.pvd 与已有索引增量合并，
所以「下一批、再跑一次」只会算新文件，且 pvd 累积成完整序列。
产物一律先写 <name>.tmp 再原子改名 —— 中途被杀不会留下半截文件被误判为已完成。
多进程并行：每个文件一个任务，-j 控制进程数，默认按可用内存自动定。

不手写任何 XML：全部交给 VTK 自己的 writer，schema 由它保证正确。

⚠️ 轴序：.fpd 里数组按 IDX = i + j*Nx + k*Nx*Ny 连续存放（x 变化最快），
   所以 numpy 必须 reshape 成 (Nz, Ny, Nx)。搞反不报错，只会得到转置的场。
   用 --self-test 验证。
"""

import argparse
import multiprocessing as mp
import os
import resource
import sys
import time
import xml.etree.ElementTree as ET

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fpd_format as F
from fpd_config import read_values


def _need_vtk():
    try:
        import vtk
        from vtk.util import numpy_support
        return vtk, numpy_support
    except ImportError:
        sys.exit("找不到 vtk 模块。本脚本需要在 pvpython 下运行：\n"
                 "  $PARAVIEW/bin/pvpython tools/fpd2vtk.py ...")


def read_arrays(path, header, want):
    """按偏移量直接读需要的数组，避免把整个文件读进内存。"""
    offs = F.field_offsets(header)
    out = {}
    with open(path, "rb") as f:
        for name in want:
            if name not in offs:
                out[name] = None
                continue
            off, n = offs[name]
            f.seek(off)
            out[name] = np.fromfile(f, dtype="<f8", count=n)
    return out


def write_atomic(path, write_fn, *args):
    """先写同目录的 <name>.tmp 再原子改名。

    直接写最终名的话，中途被打断会留下半截 .vti/.vtp，而下一轮运行只看
    「文件在不在」—— 那个半截文件就会被当成已转换而永远跳过，ParaView 里
    表现为一个静默的坏文件。改名保证「存在」== 「完整」。
    """
    tmp = path + ".tmp"
    write_fn(tmp, *args)
    os.replace(tmp, path)


def expected_products(stem, particles_only):
    """该 .fpd 对应的产物名（不含目录）。齐全即认为这个文件已转换。"""
    if particles_only:
        return [stem + ".vtp"]
    return [stem + ".vtp", stem + ".vti", stem + ".vtm"]


def mem_available():
    """MemAvailable，取自 /proc/meminfo。

    不用 os.sysconf('SC_AVPHYS_PAGES')：那个只算 free，把 buff/cache 当成已占用，
    本机（31 GB，22 GB 在 buff/cache）会报 0.6 GiB，据此定并行度只会给 1 个进程。
    """
    with open("/proc/meminfo") as f:
        for line in f:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    return None


# 单 worker 峰值内存模型：只用来定默认并行度，并在跑完后与实测 RSS 对照。
#
# 逐分量插值 + numpy_to_vtk(deep=0) 之后峰值 ≈ 36 B/格点：
#   12 B/格点   vel float32 输出
#   24 B/格点   当前那一个 float64 源场 + np.roll 副本与求和临时量
# 再加解释器 + numpy + vtk 的固定开销。取 50（比 36 多约四成余量）—— 估算偏大
# 只会少开几个进程，偏小才会把机器压垮，所以取保守侧。
#
# ⚠️ 这个常量与 write_vti 的内存写法【绑定】：优化前的写法真值是 72 B/格点
#    （三个源场、cx/cy/cz、深拷的 vel 同时活着，512×512×128 实测 2.73 GB）。
#    只改并行度不改 write_vti 的话，50 < 72 会把 jobs 抬高 —— 正是事故方向。
#
# 已知边界（不为此加参数）：忽略粒子数 N（N=28611 时仅 2.7 MB，N≈1e6 时漏
# 约 200 MB）；对 --particles-only 高估（那条路径不碰 402 MB 的场），方向安全。
MEM_BASE     = 400 << 20
MEM_PER_CELL = 50
MEM_RESERVE  = 4 << 30     # 留给桌面 —— 别再把用户的编辑器挤死


def mem_estimate(ncells):
    """单 worker 峰值内存估算（字节）。"""
    return MEM_BASE + MEM_PER_CELL * ncells


def gib(nbytes):
    """字节 -> 可读 GiB；读不到给「未知」。"""
    return "未知" if nbytes is None else "%.2f GiB" % (nbytes / 2 ** 30)


def mmss(sec):
    """秒 -> m:ss（超过一小时给 h:mm:ss）。"""
    sec = int(sec)
    if sec >= 3600:
        return "%d:%02d:%02d" % (sec // 3600, sec % 3600 // 60, sec % 60)
    return "%d:%02d" % (sec // 60, sec % 60)


def progress_line(done, total, n_new, n_skip, t0, width=24):
    """进度条单行文本（不含 \\r 等控制字符）。

    ETA 按「已完成 / 已用时间」线性外推。跳过是瞬时的、转换要几十秒，两者混在
    一起使它只是粗略值 —— 剩余文件里有多少是要真转的，读 header 之前无法知道。
    """
    el = time.time() - t0
    filled = int(width * done / total) if total else width
    eta = el * (total - done) / done if done else 0.0
    return "[%s%s] %d/%d %d%%  转换 %d  跳过 %d  已用 %s  剩余 %s" % (
        "#" * filled, "-" * (width - filled), done, total, 100 * done // total,
        n_new, n_skip, mmss(el), mmss(eta))


def default_jobs(ncells):
    """默认进程数：可用内存减掉保留量后【只用一半】，再受核数封顶。

    留一半给页缓存：输入 .fpd 是 808 MB/个、下载写进缓存转换再读出来，输出
    也有几百 MB 脏页；页缓存被挤掉时吞吐掉得比多开进程赚回来的多（实测外置盘
    热缓存下 -j 4 相对 -j 1 有 2.9×，-j 8 反而比 -j 1 还慢）。
    显式 -j 时不走这里。
    """
    avail = mem_available()
    if avail is None:
        return 1
    budget = max(0, avail - MEM_RESERVE) // 2
    return max(1, min(os.cpu_count() or 1, budget // mem_estimate(ncells)))


def convert_file(job):
    """worker：转换一个 .fpd。返回 (pvd 里的产物名, step, 是否跳过, 源文件名)。

    每个 worker 只碰自己那个 stem 的产物和 <name>.tmp，互不重名，不需要加锁。
    """
    path, out_dir, particles_only = job
    bname = os.path.basename(path)
    header = F.read_header(path)
    step = header["step"]
    stem = os.path.splitext(bname)[0]

    # 产物齐全 ⇒ 已经转过，不再读那几百 MB 重算。header（40 字节）已经读了，
    # 条目照样补进 pvd —— 上一轮若在写 pvd 之前被杀，这一轮能自愈。
    names = expected_products(stem, particles_only)
    if all(os.path.exists(os.path.join(out_dir, n)) for n in names):
        return (names[-1], step, True, bname)

    want = ["Rx", "Ry", "Rz", "Rux", "Ruy", "Ruz",
            "Vx", "Vy", "Vz", "Fx", "Fy", "Fz"]
    if not particles_only:
        want += ["vx", "vy", "vz"]
    arrs = read_arrays(path, header, want)

    vtp_name = stem + ".vtp"
    write_atomic(os.path.join(out_dir, vtp_name), write_vtp, header, arrs)
    if particles_only:
        return (vtp_name, step, False, bname)

    vti_name = stem + ".vti"
    vtm_name = stem + ".vtm"
    write_atomic(os.path.join(out_dir, vti_name), write_vti, header, arrs)
    write_atomic(os.path.join(out_dir, vtm_name), write_vtm, vti_name, vtp_name)
    return (vtm_name, step, False, bname)


def center_component(v, axis, Nz, Ny, Nx):
    """
    一个 MAC 面心速度分量插值到格心。

    v[ijk] 位于格子 ijk 在该方向上的 + 面，所以格心值 = 0.5*(v[i-1] + v[i])，
    周期边界用 np.roll 回绕。这正是 Stokes.cpp:47 已在用的约定
    （vax = 0.5*(vx[ijk] + vx[im,j,k])），直接复用，别另起一套。

    reshape 成 (Nz, Ny, Nx)：x 是最后一维（IDX 里 x 变化最快），所以
    vx/vy/vz 分别对应 axis=2/1/0。

    逐分量算而不是三个一起算，是为了让调用方算完一个就能丢掉对应的 float64
    源场：峰值只留一个 8 B/格点的输入，而不是三个。
    """
    a = v.reshape(Nz, Ny, Nx)
    return 0.5 * (a + np.roll(a, 1, axis=axis))


def stagger_to_center(vx, vy, vz, Nx, Ny, Nz):
    """三个分量一起插值。物理约定见 center_component（唯一实现）。"""
    return (center_component(vx, 2, Nz, Ny, Nx),   # x 是 axis=2
            center_component(vy, 1, Nz, Ny, Nx),   # y 是 axis=1
            center_component(vz, 0, Nz, Ny, Nx))   # z 是 axis=0


def write_vti(path, header, arrs):
    """写流体 .vti。

    ⚠️ 会 pop 掉 arrs 里的 vx/vy/vz —— 这是内存峰值所需（一次只让一个 float64
       源场活着），调用方在这之后不得再用流体场。当前唯一调用方 convert_file
       先写 .vtp（不碰这三个场）再写 .vti，顺序安全。
    """
    vtk, ns = _need_vtk()
    Nx, Ny, Nz = header["Nx"], header["Ny"], header["Nz"]

    img = vtk.vtkImageData()
    # CellData 时 WholeExtent/Dimensions 用【格点数】= 格子数 + 1
    img.SetDimensions(Nx + 1, Ny + 1, Nz + 1)
    img.SetSpacing(1.0, 1.0, 1.0)
    img.SetOrigin(0.0, 0.0, 0.0)

    # 逐分量插值，算完一个就 pop 掉它的源场：峰值只留 vel 加一个源场及其临时量。
    # 一次算出 cx/cy/cz 的写法会让三个源场和三个结果同时活着（72 B/格点）。
    vel = np.empty((Nx * Ny * Nz, 3), dtype=np.float32)
    vel[:, 0] = center_component(arrs.pop("vx"), 2, Nz, Ny, Nx).ravel()
    vel[:, 1] = center_component(arrs.pop("vy"), 1, Nz, Ny, Nx).ravel()
    vel[:, 2] = center_component(arrs.pop("vz"), 0, Nz, Ny, Nx).ravel()

    # deep=0 不深拷，省一份 12 B/格点。不悬空的依据是 numpy_support 把 numpy
    # 引用挂在了 vtkBuffer 上（GetBuffer()._numpy_reference = z），不是「vel 是
    # 局部变量」—— 所以反面是：这块内存要等 img/va 销毁才还，别把 img 提出本函数。
    va = ns.numpy_to_vtk(vel, deep=0)
    va.SetName("velocity")
    img.GetCellData().AddArray(va)
    img.GetCellData().SetVectors(va)

    w = vtk.vtkXMLImageDataWriter()
    w.SetFileName(path)
    w.SetInputData(img)
    w.SetDataModeToAppended()
    w.EncodeAppendedDataOff()
    w.Write()


def write_vtp(path, header, arrs):
    vtk, ns = _need_vtk()
    N = header["N"]

    pts = vtk.vtkPoints()
    pd = vtk.vtkPolyData()

    if N > 0:
        xyz = np.empty((N, 3), dtype=np.float64)
        xyz[:, 0] = arrs["Rx"]; xyz[:, 1] = arrs["Ry"]; xyz[:, 2] = arrs["Rz"]
        pts.SetData(ns.numpy_to_vtk(xyz, deep=1))
    pd.SetPoints(pts)

    verts = vtk.vtkCellArray()
    for i in range(N):
        verts.InsertNextCell(1)
        verts.InsertCellPoint(i)
    pd.SetVerts(verts)

    def add(name, kx, ky, kz):
        v = np.empty((N, 3), dtype=np.float32)
        v[:, 0] = arrs[kx]; v[:, 1] = arrs[ky]; v[:, 2] = arrs[kz]
        a = ns.numpy_to_vtk(v, deep=1)
        a.SetName(name)
        pd.GetPointData().AddArray(a)
        return a

    if N > 0:
        vel = add("velocity", "Vx", "Vy", "Vz")
        add("force", "Fx", "Fy", "Fz")
        # 不折叠位置：粒子穿越周期边界时 R 会瞬移，Ru 不会，做轨迹/MSD 用这个
        add("unwrapped", "Rux", "Ruy", "Ruz")
        pd.GetPointData().SetVectors(vel)

    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName(path)
    w.SetInputData(pd)
    w.SetDataModeToAppended()
    w.EncodeAppendedDataOff()
    w.Write()


def write_vtm(path, vti_name, vtp_name):
    """多块：ParaView 打开它就同时看到流体和粒子。"""
    with open(path, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="vtkMultiBlockDataSet" version="1.0" '
                'byte_order="LittleEndian" header_type="UInt32">\n')
        f.write('  <vtkMultiBlockDataSet>\n')
        f.write('    <DataSet index="0" name="fluid" file="%s"/>\n' % vti_name)
        f.write('    <DataSet index="1" name="particles" file="%s"/>\n' % vtp_name)
        f.write('  </vtkMultiBlockDataSet>\n')
        f.write('</VTKFile>\n')


def write_pvd(path, entries):
    """时间序列。timestep 用物理时间 step*dt。原子写，理由同 write_atomic。"""
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n')
        f.write('  <Collection>\n')
        for t, name in entries:
            f.write('    <DataSet timestep="%.17g" group="" part="0" file="%s"/>\n' % (t, name))
        f.write('  </Collection>\n')
        f.write('</VTKFile>\n')
    os.replace(tmp, path)


def read_pvd(path):
    """读回已有 .pvd 的 {产物名: 物理时间}。不存在返回空字典（首次运行）。"""
    if not os.path.exists(path):
        return {}
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as e:
        sys.exit("%s 不是合法 XML（%s）。\n删掉它即可重建索引，已转好的产物不受影响。"
                 % (path, e))
    return {ds.get("file"): float(ds.get("timestep")) for ds in root.iter("DataSet")}


def self_test():
    """
    判据 9/10：轴序与插值方向。
    构造 vx = 已知线性场，检查插值后的格心值与解析值一致。
    """
    Nx, Ny, Nz = 8, 4, 2
    # v[IDX(i,j,k)] = i*1e6 + j*1e3 + k
    lin = np.zeros(Nx * Ny * Nz)
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                lin[i + j * Nx + k * Nx * Ny] = i * 1e6 + j * 1e3 + k

    a = lin.reshape(Nz, Ny, Nx)
    ok_axis = (a[1, 2, 3] == 3e6 + 2e3 + 1)
    print("  轴序 reshape(Nz,Ny,Nx): a[k=1,j=2,i=3] = %.0f  期望 %.0f  %s"
          % (a[1, 2, 3], 3e6 + 2e3 + 1, "PASS" if ok_axis else "FAIL"))

    # vx 只随 i 线性变化，格心值应为 i*1e6 - 0.5e6（周期边界处除外）
    vx = np.zeros(Nx * Ny * Nz)
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                vx[i + j * Nx + k * Nx * Ny] = i * 1.0
    cx, _, _ = stagger_to_center(vx, vx, vx, Nx, Ny, Nz)
    # 内部格点 i>=1: 0.5*(i + (i-1)) = i - 0.5
    inner = cx[0, 0, 1:]
    want = np.arange(1, Nx) - 0.5
    ok_interp = np.allclose(inner, want)
    print("  x 插值 内部格点: %s  期望 %s  %s"
          % (inner, want, "PASS" if ok_interp else "FAIL"))

    return 0 if (ok_axis and ok_interp) else 1


def main():
    ap = argparse.ArgumentParser(description=".fpd -> VTK（离线转换）")
    ap.add_argument("files", nargs="*", help=".fpd 文件，可用通配符")
    ap.add_argument("-o", "--out", default="vis", help="输出目录")
    ap.add_argument("--stride", type=int, default=1, help="每隔几个文件取一个")
    ap.add_argument("-j", "--jobs", type=int, default=0,
                    help="并行进程数（0 = 按可用内存与核数自动定）")
    ap.add_argument("--config", help="提供 dt 的模拟配置文件")
    ap.add_argument("--particles-only", action="store_true")
    ap.add_argument("--self-test", action="store_true", help="轴序与插值方向自检")
    a = ap.parse_args()

    if a.self_test:
        return self_test()

    if not a.files:
        ap.error("需要至少一个 .fpd 文件")

    if not a.config:
        ap.error("需要 --config")
    dt = read_values(a.config, ("dt",))["dt"]
    files = sorted(a.files)[:: a.stride]
    if not os.path.isdir(a.out):
        os.makedirs(a.out)

    # base_name 与并行度都用第一个文件推：同一批的网格与 N 相同。
    stem0 = os.path.splitext(os.path.basename(files[0]))[0]
    base_name = stem0.rsplit("_", 1)[0] if "_" in stem0 else stem0
    ncells = F.size_of(F.read_header(files[0]))
    jobs = a.jobs if a.jobs > 0 else default_jobs(ncells)
    jobs = min(jobs, len(files))        # 别开比任务还多的进程
    print("并行 %d 个进程，共 %d 个文件（可用内存 %s，单文件估算 %s）"
          % (jobs, len(files), gib(mem_available()), gib(mem_estimate(ncells))))

    entries = {}          # {产物名: 物理时间}，本次转换的和跳过的都记进来
    n_new = n_skip = done = 0
    # 交互终端上画单行进度条原地刷新；重定向到日志时改回逐文件一行 —— 日志里塞
    # \r 控制字符会把日志毁掉，而这些行本来就是日志要看的东西。
    live = sys.stdout.isatty()
    t0 = time.time()
    work = [(p, a.out, a.particles_only) for p in files]
    # 用 spawn 而不是默认的 fork：numpy 的 OpenBLAS 在 import 时就建了线程池，
    # 父进程因而是多线程的，而 fork 只复制调用线程 —— Python 3.12 会直接报
    # 「use of fork() may lead to deadlocks in the child」。本仓库一贯不赌这类
    # 无诊断故障。代价是每个 worker 重新 import 一次 numpy+vtk（各约 0.3 s），
    # 对 808 MB/个 的文件可忽略；换来不继承父进程的线程与 VTK 全局状态。
    #
    # 已知残余风险（不加超时层）：worker 被 OOM killer 杀掉时，imap_unordered
    # 会等一个永远不来的结果 —— Pool 会补 worker 但不会补任务，表现为【挂住】
    # 而不是报错。所以上面那行内存预估该看。
    ctx = mp.get_context("spawn")
    with ctx.Pool(jobs) as pool:
        # imap_unordered：谁先转完谁先报，不等前面的慢文件。chunksize=1 是因为
        # 单个文件本来就够大，再批就是纯拖延。
        for name, step, skipped, bname in pool.imap_unordered(convert_file, work, chunksize=1):
            entries[name] = step * dt
            done += 1
            if skipped:
                n_skip += 1
            else:
                n_new += 1
            if live:
                # 行尾 \x1b[K 擦到行末：秒数从 0:05 长到 1:00:00 或缩短时不留残影。
                sys.stdout.write("\r" + progress_line(done, len(files), n_new, n_skip, t0)
                                 + "\x1b[K")
            elif skipped:
                sys.stdout.write("  跳过 %s（产物已存在）\n" % bname)
            else:
                sys.stdout.write("  %s -> %s%s\n"
                                 % (bname, name, "" if a.particles_only else " (+ .vti/.vtp)"))
            sys.stdout.flush()
    if live and done:
        sys.stdout.write("\n")      # 给进度条那行收尾，免得后续输出挤在同一行

    # 单 worker 峰值 RSS 实测。Linux 的 ru_maxrss 单位是 KiB，取所有已回收子进程
    # 里最大的那个；必须在 with 块之后读，Pool 退出时会 join 全部 worker。
    # 与估算并列，让内存模型每次运行都自校 —— 上次事故之所以无声，就是因为
    # 没有任何地方报告内存。n_new==0 时不报：那种情况 worker 只读了 40 字节的
    # header、连 vtk 都没 import，报出来会把校准带偏。
    if n_new:
        peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * 1024
        print("单 worker 峰值 RSS 实测 %s（模型估算 %s）" % (gib(peak), gib(mem_estimate(ncells))))

    # 与已有 .pvd 合并：本次没扫到的条目（别的批次、被 --stride 跳过的）保留。
    pvd = os.path.join(a.out, base_name + ".pvd")
    merged = read_pvd(pvd)
    merged.update(entries)
    ordered = sorted(merged.items(), key=lambda kv: (kv[1], kv[0]))
    write_pvd(pvd, [(t, name) for name, t in ordered])
    print("\n扫描 %d 个文件：转换 %d，跳过 %d。pvd 现含 %d 个时间步。"
          "在 ParaView 里打开：\n  %s" % (len(files), n_new, n_skip, len(ordered), pvd))
    return 0


if __name__ == "__main__":
    sys.exit(main())
