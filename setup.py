from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='diffusion_kinetics',
      version='1.0',
      description='A diffusion kinetics code for simulating the chemistry of planetary atmospheres',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
      ],
      keywords=['planet','atmospheres'],
#      url='http://github.com/storborg/funniest',
      author='Dong Wang',
      author_email='wangdong8500@gmail.com',
      license='MIT',
      packages=['diffusion_kinetics'],
      install_requires=[
          'numpy'
      ],
      include_package_data=True,
      zip_safe=False
      )


