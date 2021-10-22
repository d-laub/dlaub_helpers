A package of helper functions in Python, currently for bioinformatics work.

### Installation
Download the repo from GitHub and add the location to your PYTHONPATH environment variable. For example:
```python
repo_path='$HOME/dlaub/.program/dlaub_helpers'
export PYTHONPATH=${repo_path}${PATH:+:${PATH}}
```
Then make sure you have the requirements installed. See `requirements.txt`.

As this package grows I will improve the installation process and documentation.

### Functions
`helpers.rnaseq.pl.volcano()`

Creates volcano plots like this one.

<img src="https://s3.us-west-2.amazonaws.com/secure.notion-static.com/f66b4079-2b1e-46e7-9a44-62067335f54e/Fig4_PanelA.png?X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Credential=AKIAT73L2G45O3KS52Y5%2F20211022%2Fus-west-2%2Fs3%2Faws4_request&X-Amz-Date=20211022T194604Z&X-Amz-Expires=86400&X-Amz-Signature=e76c46b60d22d7094188b5f3ffb67b2818875db1d9ae9db34e547763e5bfafba&X-Amz-SignedHeaders=host&response-content-disposition=filename%20%3D%22Fig4_PanelA.png%22" width="300"/>
