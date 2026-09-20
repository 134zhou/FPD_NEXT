"""离线工具读取模拟配置中的数值参数；纯 stdlib。"""


def read_values(path, keys):
    values = {}
    with open(path, encoding="utf-8") as f:
        for lineno, line in enumerate(f, 1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                raise ValueError("%s:%d 不是 key=value" % (path, lineno))
            key, value = (part.strip() for part in line.split("=", 1))
            if key in keys:
                values[key] = float(value)
    missing = set(keys) - values.keys()
    if missing:
        raise ValueError("%s 缺少 %s" % (path, ", ".join(sorted(missing))))
    return values
