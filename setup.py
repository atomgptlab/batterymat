import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="batterymat",
    version="2024.2.5",
    author="Kamal Choudhary, Michael Woodcox",
    author_email="kamal.choudhary@nist.gov",
    description="batterymat",
    install_requires=[
        "numpy>=1.22.0",
        "scipy>=1.6.3",
        "jarvis-tools>=2021.07.19",
        "alignn",
        "bokeh",
        "bokeh_sampledata",
        "matplotlib",
    ],
    #scripts=["chemnlp/utils/run_chemnlp.py"],
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/atomgptlab/batterymat",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
