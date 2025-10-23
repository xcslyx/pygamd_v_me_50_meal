from setuptools import setup, find_packages

setup(
    name="pygamd_v_me_50_meal",
    version="0.1.0",
    packages=find_packages(),
    install_requires=["numpy", "pandas", "matplotlib",
                        "torch", "tqdm", "scipy",
                        "scikit-learn", "MDAnalysis", ],  # 如果有依赖包可以在这里写，例如 ["numpy", "pandas"]
    entry_points={
        'console_scripts': [
            'v50 = __main__:main'
        ],
    },
    python_requires=">=3.11",
)
