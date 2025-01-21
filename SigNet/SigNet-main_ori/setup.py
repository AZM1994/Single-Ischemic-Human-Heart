from setuptools import setup, find_packages

# Read the content of your README file
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='signaturesnet',
      version='0.1.1',
      description="Package to manipulate mutational processes.",
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/weghornlab/SigNet',
      packages=find_packages(),
      install_requires=[
            'torch',
            'scipy',
            'numpy',
            'matplotlib',
            'pandas',
            'seaborn',
            'scikit_optimize',
            'tqdm',
            'pyparsing',
            # 'gaussian_process',
            'PyYAML',
            'scikit_learn',
            'openpyxl',
            'tensorboard',
            'wandb',
      ],
      package_data={
            'signaturesnet': [
                  'data/data.xlsx',
                  'data/mutation_type_order.xlsx',
                  'data/real_data/*.txt'
                  'data/real_data/*.csv'
                  'trained_models/**/*.json',
                  'trained_models/**/*.zip',
                  'configs/**/*.yaml'
            ],
      },
)
