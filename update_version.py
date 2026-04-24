#!/usr/bin/env python3
"""
更新项目版本号的脚本
该脚本会同时更新pyproject.toml和__init__.py中的版本号
"""
import os
import re
import sys

def update_version(new_version):
    """更新版本号"""
    # 更新pyproject.toml
    pyproject_path = os.path.join(os.path.dirname(__file__), 'pyproject.toml')
    with open(pyproject_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 替换版本号
    new_content = re.sub(r'version = "[0-9]+\.[0-9]+\.[0-9]+"', f'version = "{new_version}"', content)
    
    with open(pyproject_path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    # 更新__init__.py
    init_path = os.path.join(os.path.dirname(__file__), 'pygamd_v_me_50_meal', '__init__.py')
    with open(init_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 替换版本号
    new_content = re.sub(r'__version__ = "[0-9]+\.[0-9]+\.[0-9]+"', f'__version__ = "{new_version}"', content)
    
    with open(init_path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    
    print(f"版本号已更新为: {new_version}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python update_version.py <新版本号>")
        print("例如: python update_version.py 0.5.5")
        sys.exit(1)
    
    new_version = sys.argv[1]
    # 验证版本号格式
    if not re.match(r'^[0-9]+\.[0-9]+\.[0-9]+$', new_version):
        print("错误: 版本号格式不正确，应为 X.Y.Z 格式")
        sys.exit(1)
    
    update_version(new_version)
