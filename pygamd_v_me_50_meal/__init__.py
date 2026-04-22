import tomli
import os

# 从pyproject.toml中读取版本号
with open(os.path.join(os.path.dirname(__file__), "..", "pyproject.toml"), "rb") as f:
    pyproject_data = tomli.load(f)
    __version__ = pyproject_data["project"]["version"]
