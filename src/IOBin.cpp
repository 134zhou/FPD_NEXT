#include "IOBin.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <vector>

// ---------------------------------------------------------------------------
// FNV-1a 64 位。Python 侧（tools/make_init.py）必须是逐位相同的实现 ——
// 两处独立实现同一个哈希是典型漂移点，改动任一侧都要重跑互操作判据。
// ---------------------------------------------------------------------------
static const uint64_t FNV_OFFSET = 14695981039346656037ULL;
static const uint64_t FNV_PRIME  = 1099511628211ULL;

struct Hasher
{
    uint64_t h = FNV_OFFSET;

    void feed(const void* buf, size_t n)
    {
        const unsigned char* p = (const unsigned char*)buf;
        for (size_t i = 0; i < n; i++)
        {
            h ^= (uint64_t)p[i];
            h *= FNV_PRIME;
        }
    }
};

// 写并顺带哈希
struct Writer
{
    FILE*  f;
    Hasher hash;
    bool   ok = true;

    void put(const void* buf, size_t n)
    {
        if (!ok) { return; }
        if (fwrite(buf, 1, n, f) != n) { ok = false; return; }
        hash.feed(buf, n);
    }
    template<class T> void val(const T& v) { put(&v, sizeof(T)); }
    void arr(const double* a, size_t n) { put(a, n * sizeof(double)); }
};

// 读并顺带哈希
struct Reader
{
    FILE*  f;
    Hasher hash;
    bool   ok = true;

    void get(void* buf, size_t n)
    {
        if (!ok) { return; }
        if (fread(buf, 1, n, f) != n) { ok = false; return; }
        hash.feed(buf, n);
    }
    template<class T> void val(T& v) { get(&v, sizeof(T)); }
    void arr(double* a, size_t n) { get(a, n * sizeof(double)); }
};

// header 的字段顺序在写和读之间必须严格一致，所以用同一个宏展开两次
#define FPD_HEADER_FIELDS(OP, h)  \
    OP(h.version);                                                \
    OP(h.Nx); OP(h.Ny); OP(h.Nz); OP(h.N);                          \
    OP(h.step);

std::string ckpt_path(const std::string& out_dir, const std::string& run_name, long step)
{
    char buf[32];
    snprintf(buf, sizeof(buf), "_%07ld.fpd", step);
    return out_dir + "/" + run_name + buf;
}

bool ensure_dir(const std::string& path, std::string& err)
{
    if (path.empty()) { return true; }
    struct stat st;
    if (stat(path.c_str(), &st) == 0)
    {
        if (S_ISDIR(st.st_mode)) { return true; }
        err = path + " 已存在但不是目录";
        return false;
    }
    if (mkdir(path.c_str(), 0755) != 0)
    {
        err = std::string("无法创建目录 ") + path;
        return false;
    }
    return true;
}

bool save_checkpoint(const char* path, const CkptHeader& h_in,
                     const CkptArrays& a, std::string& err)
{
    CkptHeader h = h_in;

    FILE* f = fopen(path, "wb");
    if (!f) { err = std::string("无法写入 ") + path; return false; }

    Writer w; w.f = f;

    char magic[8]; memset(magic, 0, 8); memcpy(magic, FPD_MAGIC, 7);
    w.put(magic, 8);
    uint32_t endian = FPD_ENDIAN_TAG;
    w.val(endian);

    #define WOP(x) w.val(x)
    FPD_HEADER_FIELDS(WOP, h)
    #undef WOP

    const size_t size = h.size();
    const size_t N    = (size_t)h.N;

    w.arr(a.vx, size); w.arr(a.vy, size); w.arr(a.vz, size);

    w.arr(a.Rx, N);  w.arr(a.Ry, N);  w.arr(a.Rz, N);
    w.arr(a.Rux, N); w.arr(a.Ruy, N); w.arr(a.Ruz, N);
    w.arr(a.Vx, N);  w.arr(a.Vy, N);  w.arr(a.Vz, N);
    w.arr(a.Fx, N);  w.arr(a.Fy, N);  w.arr(a.Fz, N);

    // 校验和本身不参与哈希
    const uint64_t sum = w.hash.h;
    if (w.ok && fwrite(&sum, 1, 8, f) != 8) { w.ok = false; }

    if (fclose(f) != 0) { w.ok = false; }
    if (!w.ok)
    {
        err = std::string("写入 ") + path + " 失败（磁盘满？）";
        return false;
    }
    return true;
}

// 读 header 并把 FILE* 停在数组区起点；ok 时返回打开的句柄
static FILE* open_and_read_header(const char* path, CkptHeader& h,
                                  Reader& r, std::string& err)
{
    FILE* f = fopen(path, "rb");
    if (!f) { err = std::string("打不开 ") + path; return 0; }

    r.f = f;

    char magic[8];
    r.get(magic, 8);
    if (!r.ok || memcmp(magic, FPD_MAGIC, 7) != 0)
    {
        err = std::string(path) + " 不是 .fpd 文件（magic 不匹配）";
        fclose(f); return 0;
    }

    uint32_t endian = 0;
    r.val(endian);
    if (!r.ok || endian != FPD_ENDIAN_TAG)
    {
        std::ostringstream o;
        o << path << " 的字节序标记是 0x" << std::hex << endian
          << "，本机期望 0x" << FPD_ENDIAN_TAG << "（跨架构的文件不支持）";
        err = o.str();
        fclose(f); return 0;
    }

    #define ROP(x) r.val(x)
    FPD_HEADER_FIELDS(ROP, h)
    #undef ROP

    if (!r.ok) { err = std::string(path) + " 头部被截断"; fclose(f); return 0; }

    if (h.version != FPD_VERSION)
    {
        std::ostringstream o;
        o << path << " 的版本是 " << h.version << "，本程序只认 " << FPD_VERSION;
        err = o.str();
        fclose(f); return 0;
    }
    if (h.Nx <= 0 || h.Ny <= 0 || h.Nz <= 0 || h.N < 0)
    {
        err = std::string(path) + " 头部的维度不合法";
        fclose(f); return 0;
    }
    return f;
}

bool open_checkpoint(const char* path, CkptHeader& h, CkptReader& ctx, std::string& err)
{
    Reader r;
    FILE* f = open_and_read_header(path, h, r, err);
    if (!f) { return false; }
    ctx.f = f;
    ctx.hash = r.hash.h;
    return true;
}

void close_checkpoint(CkptReader& ctx)
{
    if (ctx.f) { fclose((FILE*)ctx.f); ctx.f = 0; }
}

bool read_ckpt_header(const char* path, CkptHeader& h, std::string& err)
{
    CkptReader ctx;
    if (!open_checkpoint(path, h, ctx, err)) { return false; }
    close_checkpoint(ctx);
    return true;
}

bool read_checkpoint(CkptReader& ctx, const CkptHeader& h,
                     const CkptArrays& a, std::string& err)
{
    FILE* f = (FILE*)ctx.f;
    Reader r; r.f = f; r.hash.h = ctx.hash;
    const size_t size = h.size();
    const size_t N    = (size_t)h.N;

    r.arr(a.vx, size); r.arr(a.vy, size); r.arr(a.vz, size);
    r.arr(a.Rx, N);  r.arr(a.Ry, N);  r.arr(a.Rz, N);
    r.arr(a.Rux, N); r.arr(a.Ruy, N); r.arr(a.Ruz, N);
    r.arr(a.Vx, N);  r.arr(a.Vy, N);  r.arr(a.Vz, N);
    r.arr(a.Fx, N);  r.arr(a.Fy, N);  r.arr(a.Fz, N);

    if (!r.ok)
    {
        err = "检查点数据区被截断";
        close_checkpoint(ctx); return false;
    }

    uint64_t want = 0;
    const bool got = (fread(&want, 1, 8, f) == 8);
    char extra;
    const bool trailing = (fread(&extra, 1, 1, f) == 1);
    close_checkpoint(ctx);

    if (!got)     { err = "检查点缺少尾部校验和"; return false; }
    if (trailing) { err = "检查点尾部有多余数据"; return false; }
    if (want != r.hash.h)
    {
        std::ostringstream o;
        o << "检查点校验和不匹配：期望 0x" << std::hex << want
          << "，实际 0x" << r.hash.h;
        err = o.str();
        return false;
    }
    return true;
}

bool load_checkpoint(const char* path, CkptHeader& h,
                     const CkptArrays& a, std::string& err)
{
    CkptReader ctx;
    if (!open_checkpoint(path, h, ctx, err)) { return false; }
    return read_checkpoint(ctx, h, a, err);
}
