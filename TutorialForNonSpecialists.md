# Details #

CovalentDock is a molecular docking package designed for covalent bond forming ligands. Comparing with other docking softwares CovalentDock is relatively easy to use. However, it may still be not so clear for user who only uses it occasionally, for example,the experts majored on organic synthesis. I'm not explaining how to install Amber or some other things here, if you have problems about that, please google and you can found tons of solutions.

After you have configured Amber and Ambertools correctly you can start configuring CovalentDock. First unzip the package to a new folder which will then become the CovalentDock home folder,aka CD\_HOME.

Before you are ready to perform a covalent docking, you should put 3 related files in to a new folder. The 3 files are respectively the PDB file for the receptor, the MOL2 file for the ligand and the center of the binding site. The PDB file should be renamed as `'<your task name>_receptor.pdb'` and the mol2 file shoulde be renamed as `'<your task name>_ligand.mol2'`. Notice that the `<your task name>` in both files should be exactly the same and there's a under-strike `'_'` after it. The center file should just be named as 'center', containing the coordinates of the binding site separated by space.

After all the previous job was done we can run the preprocessing. Enter the folder containing the task files and type in the command:

$CD\_HOME/preprocessing '<your task name>'

and you should see the preprocessing running. When it is finished type in

cd bond

to enter the folder containing generated gpf and dpf files.Then we should run CovalentGrid just like Autogrid:

$CD\_HOME/CovalentGrid -p '<your grid parameter file>.gpf' -l '<your grid log file>.glg'

It may take a much longer time compared with Autodock. After that we could finally run CovalentDock:

$CD\_HOME/CovalentDock -p '<your docking parameter file>.dpf' -l '<your docking log file>.dlg'

You should also be careful that there may be more than one dpf files which need to run CovalentDock, and the number of generated dpf files depends on the number of new chiral centers after bonding.

After you have obtained the .dlg files you can analyze them by the way you are familiar with. If you still have any problem you can contact us via email or gtalk.

Enjoy!