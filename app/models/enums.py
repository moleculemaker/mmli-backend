from enum import Enum


class JobType(str, Enum):
    CLEAN = 'clean'
    CHEMSCRAPER = 'chemscraper'
    MOLLI = 'molli'
    NOVOSTOIC_OPTSTOIC= 'novostoic-optstoic'
    NOVOSTOIC_NOVOSTOIC = 'novostoic-novostoic'
    NOVOSTOIC_ENZRANK = 'novostoic-enzrank'
    NOVOSTOIC_DGPREDICTOR = 'novostoic-dgpredictor'
    SOMN = 'somn'
    DEFAULT = 'defaults'

    @staticmethod
    def values():
        # provide a meaningful named getter
        return JobType._value2member_map_

    def __str__(self) -> str:
        return self.value


class JobStatus(str, Enum):
    QUEUED = 'queued'
    PROCESSING = 'processing'
    COMPLETED = 'completed'
    ERROR = 'error'

    def __str__(self) -> str:
        return self.value


JobTypes = [str(job_type) for job_type in JobType]
