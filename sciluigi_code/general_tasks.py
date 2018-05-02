import os
import sciluigi as sl

class LoadFile(sl.ExternalTask):
    path = sl.Parameter()

    def out_file(self):
        return sl.ContainerTargetInfo(self, self.path)
