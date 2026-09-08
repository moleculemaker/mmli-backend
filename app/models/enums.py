from enum import Enum


class JobType(str, Enum):
    ACERETRO = 'aceretro'
    CLEAN = 'clean'
    CLEANDB_MEPESM = 'cleandb-mepesm'
    CHEMSCRAPER = 'chemscraper'
    CRISPR_COPIES = 'crispr-copies'
    MOLLI = 'molli'
    MUTAGENESIS = 'mutagenesis'
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
    EZSPEC_UNIDOCK = 'ezspec-unidock'
    EZSPEC_INFERENCE = 'ezspec-inference'
    ML_SIMPLEFOLD = 'ml-simplefold'

    def __str__(self) -> str:
        return self.value


class JobStatus(str, Enum):
    QUEUED = 'queued'
    PROCESSING = 'processing'
    COMPLETED = 'completed'
    ERROR = 'error'
    # Kubernetes has no concept of a canceled Job - a deleted one simply stops
    # existing - so this phase is only ever set by us, and the watcher must not
    # derive over the top of it.
    CANCELED = 'canceled'

    def __str__(self) -> str:
        return self.value


JobTypes = [str(job_type) for job_type in JobType]
ExampleJobTypes = [JobType.DEFAULT]
