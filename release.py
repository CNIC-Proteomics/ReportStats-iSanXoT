# import modules
import os
import shutil
import sys

# main function
if __name__ == "__main__":
    version = sys.argv[1]
    versionName = f"ReportStats-iSanXoT-{version}"
    home = os.getcwd()
    version_dir = os.path.join(home, versionName)

    os.mkdir(version_dir)

    shutil.copytree(os.path.join(home, "ChromiumPortable"), os.path.join(version_dir, "ChromiumPortable"))
    shutil.copytree(os.path.join(home, "R-Portable"), os.path.join(version_dir, "R-Portable"))
    shutil.copytree(os.path.join(home, "ShinyApp"), os.path.join(version_dir, "ShinyApp"))
    shutil.copytree(os.path.join(home, "Launcher"), os.path.join(version_dir, "Launcher"))
    shutil.copyfile(os.path.join(home, "iSanXoT-LIMMA.bat"), os.path.join(version_dir, "iSanXoT-LIMMA.bat"))

    shutil.make_archive(versionName, "zip", home, versionName)
    shutil.rmtree(version_dir)