from setuptools import setup, find_packages
import os
import sys

# 从__init__.py中读取版本号
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'pygamd_v_me_50_meal'))
from pygamd_v_me_50_meal import __version__

try:
    setup(
        version=__version__
    )
except Exception:
    # 如果导入失败，使用默认版本号
    setup(
        version="0.5.4"
    )

