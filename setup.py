import setuptools 

setuptools.setup(     
    name="cp2kRTPtools",
    version="1.0",     
    author="Guillaume Le Breton",     
    author_email="guillaume_le_breton@live.fr",     
    description="Analysis tools for CP2K, especially for Real Time Propagation.",     
    url="https://github.com/glb96/cp2kRTPtools",     
    zip_safe=False,     
    packages=setuptools.find_packages(),         
    classifiers=[      
                 "Programming Language :: Python :: 3", 
                 "License :: OSI Approved :: GPL3", 
                 "Operating System :: OS Independent",     ],     
    python_requires='>=3.6',   
    install_requires=['numpy'],
    )
