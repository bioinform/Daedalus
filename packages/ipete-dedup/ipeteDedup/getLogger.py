import logging

def getLogger(basename=None, debug=False, **kwargs):
    """
    Set up logging in a reusable manner, usable by either a wrapper script, which will pass the logger on to childern
    (i.e. submodules) or by individual submodules being used in isolation.
    When debug mode is on, all messages go to both the log file and to console.
    When debug mode is off, the INFO, WARNING, ERROR messages are sent to file, while console only sees
    WARNING and ERROR messages.
    :param basename: str representing the path+prefix of the log file.
    :param jsonfile:
    :param debug: boolean
    :return: logger instance
    """
    if basename is None:
        logging.warning('No basename or json file supplied; log will not be written to file.')
    logging.captureWarnings(True)
    logger = logging.getLogger()
    logging_filehandler = logging.FileHandler(filename=basename + '.log')
    logging_consolehandler = logging.StreamHandler()
    if debug:
        logger.setLevel(logging.DEBUG)
        logging_filehandler.setLevel(logging.DEBUG)
        logging_consolehandler.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        logging_filehandler.setLevel(logging.INFO)
        logging_consolehandler.setLevel(logging.INFO)
    file_formatter = logging.Formatter('%(asctime)s__%(module)s__%(levelname)s__%(message)s')
    console_formatter = logging.Formatter('%(message)s')
    logging_filehandler.setFormatter(file_formatter)
    logging_consolehandler.setFormatter(console_formatter)
    if basename is not None:
        logger.addHandler(logging_filehandler)
    logger.addHandler(logging_consolehandler)
    return logger
