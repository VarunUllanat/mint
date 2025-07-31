from setuptools import find_packages, setup

with open("README.md") as f:
    readme = f.read()

extras = {
    "mint": [
        "biopython",
        "deepspeed==0.5.9",
        "dm-tree",
        "lightning==1.9.5",
        "omegaconf",
        "einops",
        "scipy",
    ]
}

setup(
    name="mint",
    version="0.1.0",
    description="Learning the language of protein-protein interactions",
    long_description=readme,
    long_description_content_type="text/markdown",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    author="Varun Ullanat",
    url="https://github.com/VarunUllanat/mint",
    license="MIT",
    extras_require=extras,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
    ],
    include_package_data=True,
    zip_safe=False,
)
