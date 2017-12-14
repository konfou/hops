__version__ = '2.0.0'


def __get_abspath__():
    import os
    return os.path.abspath(os.path.dirname(__file__))


def run_app():
    import os
    os.system("python {0}".format(__get_abspath__()))