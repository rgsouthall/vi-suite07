class Language:
    """
    This class is a simple enumeration container for the different languages implemented
    """
    English = 'en'
    Spanish = 'sp'
    French = 'fr'


# This variable is a global parameter to hold the language state for the running program
CurrentLanguage = Language.English

EnglishDictionary: dict[str, str] = {
    'ABOUT_DIALOG': 'This program was created by NREL for the United States Department of Energy.',
    'About...': 'About...',
    'All transitions completed successfully - Open run directory for transitioned file':
        'All transitions completed successfully - Open run directory for transitioned file',
    'Attempting to cancel simulation ...': 'Attempting to cancel simulation ...',
    'Cancel Run': 'Cancel Run',
    'Cannot find a matching transition tool for this idf version':
        'Cannot find a matching transition tool for this idf version',
    'Choose E+ Folder...': 'Choose E+ Folder...',
    'Choose EnergyPlus Install Root': 'Choose EnergyPlus Install Root',
    'Choose File to Update...': 'Choose File to Update...',
    'Close': 'Close',
    'Completed Transition': 'Completed Transition',
    'Copy output files to original IDF directory?': 'Copy output files to original IDF directory?',
    'Could not open run directory': 'Could not open run directory',
    'EnergyPlus Installation': 'EnergyPlus Installation',
    'EnergyPlus Version': 'EnergyPlus Version',
    'Exit': 'Exit',
    'Failed Transition': 'Failed Transition',
    'File Details: ': 'File Details: ',
    'File Path': 'File Path',
    'IDF File doesn\'t exist at path given; cannot transition':
        'IDF File doesn\'t exist or invalid E+ install; cannot transition',
    'IDF File exists, ready to go': 'IDF File exists, ready to go',
    'IDF Selection': 'IDF or LST Selection',
    'Install Details: ': 'Install Details: ',
    'Invalid Version': 'Invalid Version',
    'Keep Intermediate Versions of Files?': 'Keep Intermediate Versions of Files?',
    'Language Confirmation': 'Language Confirmation',
    'List File Version': 'List file may contain multiple versions of files',
    'Menu': 'Menu',
    'Old Version': 'Original IDF Version',
    'Open Directory': 'Open Directory',
    'Open File for Transition': 'Open File for Transition',
    'OUTPUT_PATH':
        'Choose where to place output files, either in the run directory '
        'with the transition binaries, or in the original IDF directory',
    'Program Initialized': 'Program Initialized',
    'Running Transition': 'Running Transition',
    'Selected Directory: ': 'Selected Directory: ',
    'Selected IDF: ': 'Selected Input File: ',
    'Transition Cancelled': 'Transition Cancelled',
    'Transition cancelled': 'Transition cancelled',
    'Update File': 'Update File',
    'You must restart the app to make the language change take effect.  Would you like to restart now?':
        'You must restart the app to make the language change take effect.  Would you like to restart now?'
}

SpanishDictionary: dict[str, str] = {
    'ABOUT_DIALOG': 'Este programa fue creado por el NREL para el Departamento de Energia de los Estados Unidos.',
    'About...': 'Acerca de...',
    'All transitions completed successfully - Open run directory for transitioned file':
        'Todas las transiciones completada con éxito - Abrir directorio de ejecución para el archivo de la transición',
    'Attempting to cancel simulation ...': 'Intentando cancelar la simulación ...',
    'Cancel Run': 'Cancelar Ejecutar',
    'Cannot find a matching transition tool for this idf version':
        'No se puede encontrar una herramienta de transición a juego para esta versión de la FID',
    'Choose E+ Folder...': 'Elige E+ Carpeta...',
    'Choose EnergyPlus Install Root': 'Elija raíz de instalación de EnergyPlus',
    'Choose File to Update...': 'Elegir archivo para actualizar ...',
    'Close': 'Cerca',
    'Completed Transition': 'Transición completado',
    'Copy output files to original IDF directory?': 'NEED TRANSLATION',
    'Could not open run directory': 'No se pudo abrir directorio de ejecución',
    'EnergyPlus Installation': 'EnergyPlus instalación',
    'EnergyPlus Version': 'Versión del EnergyPlus',
    'Exit': 'Salida',
    'Failed Transition': 'La transición fallida',
    'File Details: ': 'Detalles del archivo: ',
    'File Path': 'Ruta de archivo',
    'IDF File doesn\'t exist at path given; cannot transition':
        'IDF El archivo no existe en la ruta dada; no puede transición',
    'IDF File exists, ready to go': 'existe IDF del archivo, listo para ir',
    'IDF Selection': 'IDF Selección',
    'Install Details: ': 'Detalles de instalación: ',
    'Invalid Version': 'Versión inválida',
    'Keep Intermediate Versions of Files?': 'Mantener versiones intermedias de Archivos?',
    'Language Confirmation': 'Confirmar idioma',
    'List File Version': 'NEED TRANSLATION',
    'Menu': 'Menú',
    'Old Version': 'Version antigua',
    'Open File for Transition': 'Abrir archivo para la Transición',
    'Open Directory': 'NEED TRANSLATION',
    'OUTPUT_PATH': 'NEED TRANSLATION',
    'Program Initialized': 'Programa Initialized',
    'Running Transition': 'Transición corriendo',
    'Selected Directory: ': 'Directorio seleccionado: ',
    'Selected IDF: ': 'Sélection de IDF: ',
    'Transition Cancelled': 'transición Cancelado',
    'Transition cancelled': 'Transición cancelada',
    'Update File': 'Actualizar archivo',
    'You must restart the app to make the language change take effect.  Would you like to restart now?':
        'Debe reiniciar la aplicacion para que el cambio de idioma tenga efecto. Le gustaria reiniciar ahora?'
}

FrenchDictionary: dict[str, str] = {
    'ABOUT_DIALOG':
        'Ce logiciel a été créé par NREL pour le Departement de l\'Energie des Etats Unis d\'Amérique (US DOE)',
    'About...': 'A propos...',
    'All transitions completed successfully - Open run directory for transitioned file':
        'Toutes les transitions effectuées avec succès - Ouvrir le répertoire du fichier mis à jour',
    'Attempting to cancel simulation ...': 'Tentative d\'annulation de la simulation ...',
    'Cancel Run': 'Annuler l\'éxécution',
    'Choose E+ Folder...': "Choisir le dossier d'E+...",
    'Choose EnergyPlus Install Root': "Choisissez la racine d'installation d'EnergyPlus",
    'Choose File to Update...': 'Choisissez un Fichier à mettre à jour ...',
    'Close': 'Fermer',
    'Completed Transition': 'Transition terminée',
    'Copy output files to original IDF directory?': 'NEED TRANSLATION',
    'Could not open run directory': 'Impossible d\'ouvrir le répertoire d\'éxecution',
    'EnergyPlus Installation': "EnergyPlus Installation",
    'EnergyPlus Version': 'Version du EnergyPlus',
    'Exit': 'Quitter',
    'Failed Transition': 'Echec de la Transition',
    'File Details: ': "Détails du fichier : ",
    'File Path': 'Chemin du fichier',
    'IDF File doesn\'t exist at path given; cannot transition':
        'IDF fichier n\'existe pas au chemin donné; transition impossible',
    'IDF File exists, ready to go': 'Le fichier IDF existe, prêt.',
    'IDF Selection': "Sélection de l'IDF",
    'Install Details: ': "Détails de l'installation : ",
    'Invalid Version': "Version invalide",
    'Keep Intermediate Versions of Files?': 'Garder les versions intermediaires des fichiers?',
    'Language Confirmation': "Confirmer la langue",
    'List File Version': 'NEED TRANSLATION',
    'Menu': 'Menu',
    'Old Version': 'Ancienne version',
    'Open File for Transition': 'Ouvrir un fichier pour la transition',
    'Open Directory': 'NEED TRANSLATION',
    'OUTPUT_PATH': 'NEED TRANSLATION',
    'Cannot find a matching transition tool for this idf version':
        'Impossible de trouver un utilitaire de Transition pour cette version d\'IDF',
    'Program Initialized': 'Programme initialisé',
    'Running Transition': 'Transition en cours',
    'Selected Directory: ': "Répertoire sélectionné : ",
    'Selected IDF: ': "IDF selectionné : ",
    'Transition Cancelled': 'Transition annulée',
    'Transition cancelled': 'Transition annulée',
    'Update File': 'Mettre à jour',
    'You must restart the app to make the language change take effect.  Would you like to restart now?':
        'Vous devez relancer le logiciel pour effectuer le changement de langue. Voulez-vous relancer maintenant?'
}


def set_language(lang):
    """
    This is the interface for changing the language, call this, save settings, then restart the program
    :param lang: A language identifier from the :py:class:`Languages` enumeration class
    """
    global CurrentLanguage
    CurrentLanguage = lang


def report_missing_keys(mute: bool = False) -> bool:
    """
    This function simply scans dictionaries to see if any keys are missing from them compared to a baseline.
    The baseline is currently the English dictionary.
    This function simply reports to the terminal.
    """
    base_keys = EnglishDictionary.keys()
    known_dictionaries = {
        'Spanish': SpanishDictionary, 'French': FrenchDictionary
    }
    any_missing = False
    for dict_name, dictionary in known_dictionaries.items():  # add here
        if not mute:  # pragma: no cover
            print("Processing missing keys from dictionary: " + dict_name)
        for key in base_keys:
            # this should never happen in unit tests, so not covering
            if key not in dictionary:  # pragma: no cover
                if not mute:
                    print("Could not find key: \"%s\"" % key)
                any_missing = True
    return True if any_missing else False


def translate(key, mute: bool = False):
    """
    This function translates a string into a dictionary.

    :param key: The string to translate
    :param mute: A unit test override to hush output messaging.
    :return: The translated string
    """
    # if for some reason blank, just return blank
    if key is None or key == "":
        return ""

    # start with English, but switch based on language
    dictionary = EnglishDictionary
    if CurrentLanguage == Language.Spanish:
        dictionary = SpanishDictionary
    elif CurrentLanguage == Language.French:
        dictionary = FrenchDictionary

    # if the key is there, return it, otherwise return a big flashy problematic statement
    if key in dictionary:
        return dictionary[key]
    else:
        if not mute:  # pragma: no cover
            print("Could not find this key in the dictionary: \"%s\"" % key)
        return "TRANSLATION MISSING"
