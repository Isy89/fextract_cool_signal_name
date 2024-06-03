from setuptools import setup, find_packages

with open("requirements.txt", "r", encoding="utf-8") as f:
    requires = []
    for line in f:
        req = line.split("#", 1)[0].strip()
        requires.append(req)

setup(
    name="fextract_cool_signal_name",
    install_requires=requires,
    entry_points={"lbfextract": [
        "cool_signal_name = fextract_cool_signal_name.plugin:hook"
    ],
        "lbfextract_cli": [
            "extract_cool_signal_name = fextract_cool_signal_name.plugin:hook_cli"
        ]
    },
    packages=find_packages('src'),
    package_dir={
        '': 'src',
    },
)