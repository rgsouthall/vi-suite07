class EPLaunchFileException(Exception):
    def __init__(self, file_path, message=''):
        super(Exception, self).__init__(self, message)
        self.file_path = file_path
