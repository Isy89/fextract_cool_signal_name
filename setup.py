from setuptools import setup, find_packages

setup(
    name="fextract_cool_signal_name",
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