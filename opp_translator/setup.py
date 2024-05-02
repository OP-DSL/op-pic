from setuptools import setup

setup(
    name="opp-translator",
    version="2.0.0",
    packages=["opp-translator"],
    install_requires=open("requirements.txt").read().splitlines(),
    python_requires=">=3.8",
    #entry_points={"console_scripts": ["opp-translator = opp-translator.__main__:main"]},
)