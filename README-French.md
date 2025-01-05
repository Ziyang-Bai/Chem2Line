<div align="center">
<p><img style="width:200px; height:200px;" src="https://github.com/Ziyang-Bai/Chem2Line/blob/main/lib/media/nctl.png?raw=true" alt=""></p>
<h2>Chem2Line</h2>
<p>Générez des formules de ligne clé pour les organiques sur votre ordinateur !</p>
<p><em>在你的电脑生成有机物的键线式！</em></p>
</div>

[Française Version](README-French.md) | [中文版](README.md) | [English Version](README-English.md)

# Chem2Line

![GitHub stars](https://img.shields.io/github/stars/Ziyang-Bai/Chem2Line)
![GitHub forks](https://img.shields.io/github/forks/Ziyang-Bai/Chem2Line)
![GitHub issues](https://img.shields.io/github/issues/Ziyang-Bai/Chem2Line)
![GitHub license](https://img.shields.io/github/license/Ziyang-Bai/Chem2Line)
[![Build CTLGUI Executable](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/build.yml/badge.svg)](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/build.yml)
[![Release CTLGUI](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/release.yml/badge.svg)](https://github.com/Ziyang-Bai/Chem2Line/actions/workflows/release.yml)

**Chem2Line** est un outil pour convertir des formules chimiques ou des représentations SMILES en images de ligne clé. Les utilisateurs peuvent entrer des formules chimiques ou des SMILES, et le programme générera les images de structure chimique correspondantes. Le programme prend en charge le changement de bases de données, l'affichage des informations sur les bases de données et la sauvegarde des images générées.

**Fonctionnalités**

- **Formule chimique en ligne clé** : Les utilisateurs peuvent entrer des formules chimiques ou des représentations SMILES, et le programme générera automatiquement les images de ligne clé correspondantes.
- **Gestion des bases de données** : Prend en charge le chargement de bases de données SMILES au format XML, le changement de bases de données et l'affichage des informations sur la base de données actuelle.
- **Sauvegarder les images** : Permet de sauvegarder les images de ligne clé générées au format PNG.
- **À propos du développeur** : Fournit des informations détaillées sur le développeur, la version du logiciel et la version du noyau.
- **Affichage de la progression** : Affiche une barre de progression lors du chargement de la base de données.
- **Support multilingue** : Prend en charge le chinois simplifié, l'anglais américain, le français et l'allemand.
- **Codes d'erreur** : Résoudre les problèmes n'a jamais été aussi simple !

**Installation et dépendances**

Cet outil nécessite les bibliothèques Python suivantes :

- tkinter : Pour créer des interfaces graphiques.
- PIL (ou Pillow) : Pour traiter et afficher des images.
- ctlcore : Bibliothèque de base personnalisée contenant des opérations liées aux structures chimiques.

```bash
pip install pillow
pip install tk
```

ctlcore est inclus dans le dépôt.

**Utilisation**

**1. Démarrer le programme**

Après avoir exécuté le programme, une interface principale s'affichera où vous pourrez entrer des formules chimiques ou des représentations SMILES et cliquer sur le bouton "Générer ligne clé" pour voir les résultats.

**2. Changer de base de données**

Cliquez sur "Base de données" > "Changer de base de données" dans la barre de menu, et sélectionnez une nouvelle base de données au format XML. Le programme chargera automatiquement la base de données et affichera la progression du chargement.

Si vous souhaitez appeler rapidement des bases de données couramment utilisées, veuillez les placer dans le répertoire `lib/db`. La base de données par défaut peut être modifiée via `config.xml` ou la fenêtre de configuration. Cependant, elle doit être dans le répertoire `lib/db`. La modification de `config.xml` n'est pas restreinte, mais elle n'est pas recommandée.

**2.1. Format de la base de données**

La base de données utilise le format XML. Le format standard est le suivant :

    ```xml
    <smiles_database>
        <!--name publisher description sont des éléments optionnels-->
        <name>Base de données standard</name>
        <publisher>Ziyang-Bai</publisher>
        <description>Ceci est une remarque</description>
        <compound>
            <!--name et cas peuvent être n'importe quels caractères, mais smiles et formula sont obligatoires.-->
            <name>Benzène</name>
            <smiles>c1ccccc1</smiles>
            <formula>C6H6</formula>
            <cas>71-43-2</cas>
        </compound>
    </smiles_database>
    ```

Vous pouvez ajouter ou supprimer des éléments selon vos besoins. Si vous avez une base de données au format sdf, vous pouvez utiliser sdf_converter.py pour la convertir au format XML. Cependant, certaines modifications sont nécessaires. Veuillez consulter le fichier sdf_converter.py pour plus de détails.

**2.2. Télécharger des bases de données existantes**

Voici quelques bases de données disponibles en téléchargement que vous pouvez utiliser facilement.
| Nom de la base de données | Lien de téléchargement |
|---------------------------|------------------------|
| Base de données par défaut | [Télécharger](https://github.com/Ziyang-Bai/Chem2Line/blob/main/lib/db/default_database.xml) |
| PubChem Carbon Organics | [Télécharger](https://github.com/Ziyang-Bai/Chem2Line/blob/main/databases/pubchem_carbon.xml) |

**2.3. Soumettre votre base de données**

Si vous avez collecté votre propre base de données et souhaitez la partager avec d'autres, vous pouvez suivre ces étapes :

1. Nommez votre fichier de base de données `your_database.xml`.

2. Créez un fork du dépôt et téléchargez votre fichier de base de données dans le répertoire `databases`.

3. Créez une Pull Request pour fusionner votre fichier de base de données dans le dépôt principal. N'oubliez pas d'inclure le nom et la description de votre base de données.

**3. Sauvegarder les images de ligne clé**

Cliquez sur "Fichier" > "Sauvegarder l'image de ligne clé" dans la barre de menu, choisissez un emplacement de sauvegarde et sauvegardez l'image de ligne clé générée (format PNG).

**4. Afficher les informations sur la base de données**

Cliquez sur "Base de données" > "À propos de la base de données" dans la barre de menu, et le programme affichera les informations sur la base de données actuelle.

**5. À propos du développeur**

Cliquez sur "À propos" > "Développeur" dans la barre de menu pour voir les informations sur le développeur et la version du logiciel.

**6. Soumettre votre fichier de langue**

Si vous souhaitez soumettre un nouveau fichier de langue pour Chem2Line, veuillez suivre ces étapes :

1. Créez un nouveau fichier XML dans le répertoire `lang`, nommé d'après le code de langue (par exemple, `fr_fr.xml`).
2. Définissez les traductions pour le texte de l'interface dans le fichier XML, par exemple :
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
3. Ajoutez le nouveau code de langue au fichier `config.xml`, par exemple :
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
4. Soumettez vos modifications et créez une Pull Request.

## Codes d'erreur

Pendant l'utilisation, vous pouvez rencontrer les codes d'erreur suivants :

- **1000** : Erreur inconnue
- **1001** : Aucun résultat trouvé
- **1002** : Échec de la génération
- **2000** : Erreur de chargement de la configuration
- **2001** : Erreur de sauvegarde de la configuration
- **3000** : Erreur de chargement de la langue

## Améliorations

1. __Expérience utilisateur__
    - Ajouter une validation des entrées pour les formules chimiques ou SMILES invalides et fournir des messages d'erreur spécifiques.
    - Prendre en charge l'historique des formules chimiques ou SMILES récemment utilisées.
2. __Structure du code__
    - Extraire la logique de l'interface utilisateur en fonctions indépendantes, en séparant davantage la logique de base et la logique de l'interface utilisateur.
    - Utiliser le multithreading pour gérer la génération d'images, évitant de bloquer l'interface principale.
3. __Extension de la base de données__
    - Ajouter la prise en charge de plus de formats de bases de données (par exemple, JSON, CSV).
    - Fournir des exemples de formules chimiques et de SMILES par défaut intégrés.
4. __Extension des fonctionnalités__
    - Ajouter une visualisation moléculaire en 3D.
    - Fournir une rotation et un zoom des images de ligne clé.
    - Prendre en charge l'exportation en tant que graphiques vectoriels (SVG).

- __Goulots d'étranglement de performance__ : Pour les fichiers de base de données XML volumineux, le chargement et l'analyse peuvent être retardés. Optimiser les méthodes de lecture des fichiers ou introduire un index.
- __Problèmes de dépendance__ : S'assurer que les utilisateurs ont installé RDKit et d'autres dépendances.
- __Problèmes de compatibilité__ : Les performances de Tkinter peuvent être limitées sur certains systèmes. Fournir une version de l'interface interactive basée sur le Web.

**Licence**

Ce programme utilise la licence Apache-2.0. Veuillez consulter le fichier LICENSE pour plus d'informations.

# Construction
## Versions pré-construites
Veuillez télécharger les versions pré-construites depuis la page des versions de Github.
## Construction rapide avec Github Action
Une fois que vous avez cloné votre dépôt Github, vous pouvez utiliser Github Action Quick Build pour construire votre programme.
Cliquez simplement sur le bouton "Actions" en haut à droite du dépôt, puis sélectionnez le workflow "Build CTLGUI Executable". Exécutez le workflow ou faites un nouveau commit, et attendez que la construction soit terminée.
Une fois la construction terminée, vous pouvez trouver le fichier exécutable généré sur la page "Actions".
## Construction manuelle
Alternativement, vous pouvez construire manuellement votre programme en utilisant les commandes suivantes :
```bash
git clone https://github.com/Ziyang-Bai/Chem2Line.git
cd Chem2Line
python -m venv myenv
myenv\Scripts\activate  # Sur Windows
# source myenv/bin/activate  # Sur macOS et Linux
pip install --upgrade pip wheel
pip install -r requirements.txt
pip install rdkit
pip install nuitka
nuitka --standalone --enable-plugin=tkinter --windows-icon-from-ico nctl.ico ctlgui.py
```
Cela construira Chem2Line en utilisant toutes les dépendances de votre environnement Python actuel.

# Rapport d'erreurs/Contact
Créez une Issue pour signaler des erreurs, ou utilisez le groupe Telegram [![Join Telegram Channel](https://img.shields.io/badge/Telegram-Chem2Line-2CA5E0?logo=telegram)](https://t.me/Chem2Line)
