<div align="center">
<p><img style="width:200px; height:200px;" src="https://github.com/Ziyang-Bai/Chem2Line/blob/main/lib/media/nctl.png?raw=true" alt=""></p>
<h2>Chem2Line</h2>
<p>Generate key line formulas for organics at your computer!</p>
<p><em>在你的电脑生成有机物的键线式！</em></p>
</div>

[English Version](README-English.md) | [中文版](README.md) | [Française Version](README-French.md)

# Chem2Line

![GitHub stars](https://img.shields.io/github/stars/Ziyang-Bai/Chem2Line)
![GitHub forks](https://img.shields.io/github/forks/Ziyang-Bai/Chem2Line)
![GitHub issues](https://img.shields.io/github/issues/Ziyang-Bai/Chem2Line)
![GitHub license](https://img.shields.io/github/license/Ziyang-Bai/Chem2Line)
[![Build CTLGUI Executable](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/build.yml/badge.svg)](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/build.yml)
[![Release CTLGUI](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/release.yml/badge.svg)](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/release.yml)

**Chem2Line** is a tool for converting chemical formulas or SMILES representations into key line images. Users can input chemical formulas or SMILES, and the program will generate the corresponding chemical structure images. The program supports changing databases, displaying database information, and saving generated images.

**Features**

- **Chemical Formula to Key Line**: Users can input chemical formulas or SMILES representations, and the program will automatically generate the corresponding key line images.
- **Database Management**: Supports loading SMILES databases in XML format, changing databases, and viewing current database information.
- **Save Images**: Allows saving generated key line images in PNG format.
- **About the Developer**: Provides detailed information about the developer, software version, and core version.
- **Progress Display**: Displays a progress bar when loading the database.
- **Multilingual Support**: Supports Simplified Chinese, American English, French, and German.
- **Error Codes**: Solving problems has never been easier!

**Installation and Dependencies**

This tool requires the following Python libraries:

- tkinter: For creating graphical user interfaces.
- PIL (or Pillow): For processing and displaying images.
- ctlcore: Custom core library containing operations related to chemical structures.

```bash
pip install pillow
pip install tk
```

ctlcore is included in the repository.

**Usage**

**1. Start the Program**

After running the program, a main interface will be displayed where you can input chemical formulas or SMILES representations and click the "Generate Key Line" button to see the results.

**2. Change Database**

Click "Database" > "Change Database" in the menu bar, and select a new XML format database. The program will automatically load the database and display the loading progress.

If you want to quickly call commonly used databases, please place them in the `lib/db` directory. The default database can be changed through `config.xml` or the configuration window. However, it must be in the `lib/db` directory. Modifying `config.xml` is not restricted, but it is not recommended.

**2.1. Database Format**

The database uses XML format. The standard format is as follows:

    ```xml
    <smiles_database>
        <!--name publisher description are optional elements-->
        <name>Standard Database</name>
        <publisher>Ziyang-Bai</publisher>
        <description>This is a remark</description>
        <compound>
            <!--name and cas can be any characters, but smiles and formula are required.-->
            <name>Benzene</name>
            <smiles>c1ccccc1</smiles>
            <formula>C6H6</formula>
            <cas>71-43-2</cas>
        </compound>
    </smiles_database>
    ```

You can add or delete elements as needed. If you have an sdf format database, you can use sdf_converter.py to convert it to XML format. However, some modifications are needed. Please refer to the sdf_converter.py file for details.

**2.2. Download Existing Databases**

Here are some databases available for download that you can use conveniently.
| Database Name | Download Link |
|---------------|---------------|
| Default Database | [Download](https://github.com/Ziyang-Bai/Chem2Line/blob/main/lib/db/default_database.xml) |
| PubChem Carbon Organics | [Download](https://github.com/Ziyang-Bai/Chem2Line/blob/main/databases/pubchem_carbon.xml) |

**2.3. Submit Your Database**

If you have collected your own database and want to share it with others, you can follow these steps:

1. Name your database file `your_database.xml`.

2. Create a fork of the repository and upload your database file to the `databases` directory.

3. Create a Pull Request to merge your database file into the main repository. Remember to include the name and description of your database.

**3. Save Key Line Images**

Click "File" > "Save Key Line Image" in the menu bar, choose a save location, and save the generated key line image (PNG format).

**4. View Database Information**

Click "Database" > "About Database" in the menu bar, and the program will display the current database information.

**5. About the Developer**

Click "About" > "Developer" in the menu bar to view information about the developer and software version.

**6. Submit Your Language File**

If you want to submit a new language file for Chem2Line, please follow these steps:

1. Create a new XML file in the `lang` directory, named after the language code (e.g., `fr_fr.xml`).
2. Define the translations for the interface text in the XML file, for example:
    ```xml
    <language>
        <input_label>Veuillez entrer une formule chimique ou SMILES :</input_label>
        <submit_button>Générer Bondline</submit_button>
        <save_image>Enregistrer l'image Bondline</save_image>
        <exit>Quitter</exit>
        <change_database>Changer de base de données</change_database>
        <database_info>À propos de la base de données</database_info>
        <developer>Développeur</developer>
        <repository>Référentiel</repository>
        <file>Fichier</file>
        <database>Base de données</database>
        <about>À propos</about>
        <language>Langue</language>
        <english>Anglais</english>
        <chinese>Chinois</chinese>
    </language>
    ```
3. Add the new language code to the `config.xml` file, for example:
    ```xml
    <config>
        <language>zh_cn</language>
        <available_languages>
            <language>zh_cn</language>
            <language>en_us</language>
            <language>fr_fr</language>
        </available_languages>
    </config>
    ```
4. Submit your changes and create a Pull Request.

## Error Codes

During use, you may encounter the following error codes:

- **1000**: Unknown error
- **1001**: No results found
- **1002**: Generation failed
- **2000**: Configuration load error
- **2001**: Configuration save error
- **3000**: Language load error

## Improvements

1. __User Experience__
    - Add input validation for invalid chemical formulas or SMILES and provide specific error messages.
    - Support history of recently used chemical formulas or SMILES.
2. __Code Structure__
    - Extract GUI logic into independent functions, further separating core logic and user interface logic.
    - Use multithreading to handle image generation, avoiding blocking the main interface.
3. __Database Expansion__
    - Add support for more database formats (e.g., JSON, CSV).
    - Provide built-in default chemical formulas and SMILES examples.
4. __Feature Expansion__
    - Add 3D molecular visualization.
    - Provide image rotation and zooming for key line images.
    - Support exporting as vector graphics (SVG).

- __Performance Bottlenecks__: For large XML database files, loading and parsing may be delayed. Optimize file reading methods or introduce indexing.
- __Dependency Issues__: Ensure users have installed RDKit and other dependencies.
- __Compatibility Issues__: Tkinter's performance may be limited on some systems. Provide a web-based interactive interface version.

**License**

This program uses the Apache-2.0 license. Please refer to the LICENSE file for more information.

# Build
## Pre-built Versions
Please download pre-built versions from the Github Releases page.
## Github Action Quick Build
Once you have cloned your Github repository, you can use Github Action Quick Build to build your program.
Simply click the "Actions" button at the top right of the repository, then select the "Build CTLGUI Executable" workflow. Run the workflow or make a new commit, and wait for the build to complete.
Once the build is complete, you can find the generated executable file on the "Actions" page.
## Manual Build
Alternatively, you can manually build your program using the following commands:
```bash
git clone https://github.com/Ziyang-Bai/Chem2Line.git
cd Chem2Line
python -m venv myenv
myenv\Scripts\activate  # On Windows
# source myenv/bin/activate  # On macOS and Linux
pip install --upgrade pip wheel
pip install -r requirements.txt
pip install rdkit
pip install nuitka
nuitka --standalone --enable-plugin=tkinter --windows-icon-from-ico=nctl.ico ctlgui.py
```
This will build Chem2Line using all dependencies in your current Python environment.

# Error Reporting/Contact
Create an Issue to report errors, or use the Telegram group [![Join Telegram Channel](https://img.shields.io/badge/Telegram-Chem2Line-2CA5E0?logo=telegram)](https://t.me/Chem2Line)
