# pygamd_v_me_50_meal/version_check.py

import importlib.metadata
import json
import os
import sys
import urllib.request
from datetime import datetime, timezone
from packaging.version import Version

OWNER = "lyxlyxlyxxx"
REPO = "pygamd_v_me_50_meal"
PKG_NAME = "pygamd_v_me_50_meal"


def _debug(msg: str):
    if os.getenv("PYGAMD_UPDATE_DEBUG"):
        print(f"[pygamd:update] {msg}", file=sys.stderr)


def _get_local_version():
    try:
        v = Version(importlib.metadata.version(PKG_NAME))
        _debug(f"Local version detected: {v}")
        return v
    except Exception as e:
        _debug(f"Failed to read local version: {e}")
        return None


def _get_repo_info(timeout):
    url = f"https://gitee.com/api/v5/repos/{OWNER}/{REPO}"
    _debug(f"Querying repo info: {url}")
    with urllib.request.urlopen(url, timeout=timeout) as resp:
        return json.load(resp)


# 【核心修改点】：不再请求 releases，改为请求 tags 接口并解析最大版本号
def _get_latest_tag_version(timeout):
    url = f"https://gitee.com/api/v5/repos/{OWNER}/{REPO}/tags"
    _debug("Querying tags info")

    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.load(resp)
    except Exception as e:
        _debug(f"Failed to fetch tags: {e}")
        return None

    if not data:
        _debug("No tags found in repository")
        return None

    versions = []
    for tag_info in data:
        tag_name = tag_info.get("name")
        if tag_name:
            try:
                versions.append(Version(tag_name.lstrip("v")))
            except Exception as e:
                _debug(f"Invalid tag format '{tag_name}': {e}")

    if not versions:
        return None

    latest_version = max(versions)
    _debug(f"Latest tag version: {latest_version}")
    return latest_version


def _parse_utc(ts: str):
    return datetime.fromisoformat(ts.replace("Z", "+00:00"))


def check_update(timeout: float = 2.0):
    if os.getenv("PYGAMD_NO_UPDATE_CHECK"):
        _debug("Update check disabled by environment variable")
        return

    try:
        local_version = _get_local_version()
        repo_info = _get_repo_info(timeout)

        # ===== 优先级 1：基于 Tag 的版本对比 =====
        try:
            latest_version = _get_latest_tag_version(timeout)
        except Exception as e:
            latest_version = None
            _debug(f"Tag query failed: {e}")

        if local_version and latest_version:
            if latest_version > local_version:
                print(
                    f"[pygamd] 新版本可用: "
                    f"{local_version} → {latest_version}\n"
                    f"https://gitee.com/{OWNER}/{REPO}"
                )
            else:
                _debug("Local version is up to date with latest tag")
            return

        _debug("Tag-based check not available, using pushed_at fallback")

        # ===== fallback：基于推送时间的对比 =====
        pushed_at = repo_info.get("pushed_at")
        if not pushed_at:
            _debug("No pushed_at field in repo info")
            return

        remote_time = _parse_utc(pushed_at)
        _debug(f"Repository last pushed at: {remote_time.isoformat()}")

        try:
            dist = importlib.metadata.distribution(PKG_NAME)
            local_time = datetime.fromtimestamp(
                dist._path.stat().st_mtime, tz=timezone.utc
            )
            _debug(f"Local install time: {local_time.isoformat()}")
        except Exception as e:
            _debug(f"Failed to determine local install time: {e}")
            return

        if remote_time > local_time:
            print(
                "[pygamd] 远程仓库有更新。\n"
                f"请运行 'pip install git+https://gitee.com/{OWNER}/{REPO}'"
            )
        else:
            _debug("Local installation is newer than repo push time")

    except Exception as e:
        _debug(f"Unexpected failure during update check: {e}")


if __name__ == "__main__":
    # 如果你想测试，可以取消下面这行的注释来查看详细日志
    # os.environ["PYGAMD_UPDATE_DEBUG"] = "1"
    check_update()