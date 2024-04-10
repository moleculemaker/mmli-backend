from enum import Enum


class JobType(str, Enum):
    CLEAN = 'clean'
    CHEMSCRAPER = 'chemscraper'
    MOLLI = 'molli'
    SOMN = 'somn'

    def __str__(self) -> str:
        return self.value


class JobStatus(str, Enum):
    QUEUED = 'queued'
    PROCESSING = 'processing'
    COMPLETED = 'completed'
    ERROR = 'error'

    def __str__(self) -> str:
        return self.value
