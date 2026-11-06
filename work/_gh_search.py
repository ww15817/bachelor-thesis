"""查 GitHub 上有没有可用的论文/开题写作辅助仓库。用完即删。"""

import io
import json
import sys
import urllib.parse
import urllib.request

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")


def search(q: str, n: int = 6):
    url = (
        "https://api.github.com/search/repositories?q="
        + urllib.parse.quote(q)
        + f"&sort=stars&order=desc&per_page={n}"
    )
    req = urllib.request.Request(url, headers={"Accept": "application/vnd.github+json",
                                              "User-Agent": "dsh-check"})
    with urllib.request.urlopen(req, timeout=25) as r:
        d = json.load(r)
    print(f"\n=== {q} ===")
    for it in d.get("items", []):
        desc = (it.get("description") or "").replace("\n", " ")[:80]
        print(f"{it['stargazers_count']:6d}  {it['full_name']:42s}  {desc}")


for q in [
    "毕业论文 latex 模板",
    "thesis latex template chinese university",
    "开题报告",
    "academic writing skill agent",
    "claude skill writing",
]:
    try:
        search(q)
    except Exception as exc:  # noqa: BLE001
        print(f"\n=== {q} ===\n查询失败: {exc}")
