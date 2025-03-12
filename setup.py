from setuptools import setup, find_packages

setuptools.setup()

setup(
    name="biokit",
    version="0.1",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "biokit-cli = biokit.cli:main"  # 入口指向 cli.py 的 main 函数
        ]
    },
    setup_requires=['pbr'],
    pbr=True
)
