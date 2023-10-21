from enum import Enum


class JobType(str, Enum):
    CLEAN = 'clean'
    CHEMSCRAPER = 'chemscraper'
    MOLLI = 'molli'

    def __str__(self) -> str:
        return self.value


class JobStatus(str, Enum):
    QUEUED = 'queued'
    PENDING = 'pending'
    PROCESSING = 'processing'
    READY = 'ready'

    def __str__(self) -> str:
        return self.value
