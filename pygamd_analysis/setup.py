from setuptools import setup, find_packages

setup(
    name='pygamd_analysis',
    version='0.0.20',
    packages=find_packages(),
    install_requires=[
        "numpy", "pandas", "matplotlib", "torch",
        "tqdm", "scipy", "scikit-learn",
        "pygamd",
    ],
    entry_points={
        'console_scripts': [
            'pygamd_analysis=pygamd_analysis.main:main',
        ],
    },
)
