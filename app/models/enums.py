from enum import Enum


class JobType(str, Enum):
    ACERETRO = 'aceretro'
    CLEAN = 'clean'
    CLEANDB_MEPESM = 'cleandb-mepesm'
    CHEMSCRAPER = 'chemscraper'
    MOLLI = 'molli'
    NOVOSTOIC_OPTSTOIC= 'novostoic-optstoic'
    NOVOSTOIC_PATHWAYS = 'novostoic-pathways'
    NOVOSTOIC_ENZRANK = 'novostoic-enzrank'
    NOVOSTOIC_DGPREDICTOR = 'novostoic-dgpredictor'
    REACTIONMINER = 'reactionminer'
    SOMN = 'somn'
    OED_CHEMINFO = 'oed-cheminfo'
    OED_DLKCAT = 'oed-dlkcat'
    OED_UNIKP = 'oed-unikp'
    OED_CATPRED = 'oed-catpred'
    DEFAULT = 'defaults'
    EZ_SPECIFICITY = 'ez-specificity'

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
