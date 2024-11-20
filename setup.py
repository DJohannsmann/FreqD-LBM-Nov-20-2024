import subprocess
import sys

packages_versions = {
            'numpy':           '1.26.4', 
            'numba':           '0.59.0', 
            'pandas':          '2.1.4', 
            'lmfit':           '1.2.2', 
            'tk':              '0.1.0', 
            'scipy':           '1.11.4', 
            'configparser':    '7.1.0', 
            'queuelib':        '1.6.2',
            'tooltip':         '1.0.0',
            'watchdog':        '6.0.0', 
            'matplotlib':      '3.9.2'
            }

def main(package_info: dict[str, str]) -> None:
    #getting pip
    try:
        subprocess.run(["python", "get-pip.py"])
    except Exception as E: 
        print(f'Problem in getting pip occured: {repr(E)} \n \
              please install pip manually (https://pip.pypa.io/en/stable/installation/#get-pip-py)')
        pass
    
    #getting packages
    for package, version in package_info.items():
        try:
            print(f"check for {package}=={version} and install it if necessary")
            subprocess.check_call([sys.executable, "-m", "pip", "install", f"{package}=={version}"])
        except subprocess.CalledProcessError: 
            print(f"problem with {package}=={version} -> install another suitable version")
            subprocess.check_call([sys.executable, "-m", "pip", "install", f"{package}"])
        
if __name__ == "__main__":
    main(package_info = packages_versions)
