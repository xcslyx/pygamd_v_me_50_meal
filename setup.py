from setuptools import setup, find_packages

setup(
    name="pygamd_v_me_50_meal",
    version="0.1.10",
    packages=find_packages(),
    install_requires=[
        "numpy", "pandas", "matplotlib",
        "torch", "tqdm", "scipy",
        "scikit-learn", "MDAnalysis"
    ],
    entry_points={
        "console_scripts": [
            "v50=pygamd_v_me_50_meal.__main__:main"
        ]
    },
    author="Yuxi Leng",
    description="MD analysis toolkit for PYGAMD Toolkit",
    python_requires=">=3.11",
)

