from abc import ABC, abstractmethod


class BaseDataLoader(ABC):

    @abstractmethod
    def load():
        raise NotImplementedError

    @abstractmethod
    def save():
        raise NotImplementedError


class BasePreloader(ABC):

    @abstractmethod
    def get_rules():
        raise NotImplementedError

    @abstractmethod
    def get_sink():
        raise NotImplementedError

    @abstractmethod
    def get_sources():
        raise NotImplementedError
