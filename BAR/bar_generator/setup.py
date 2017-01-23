from setuptools import setup

setup(name='bar_generator',
      version='0.1',
      description='BAR simulation automation script.',
      url='http://github.com/Nevensky/Free-Energy-Skripte',
      author='Neven Golenic',
      author_email='neven.golenic@gmail.com',
      license='MIT',
      packages=['funniest'],
      install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'termcolor',
          'click'
      ],
      zip_safe=False)