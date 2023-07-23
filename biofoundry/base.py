from abc import ABC, abstractmethod


class BaseDataLoader(ABC):

    @abstractmethod
    def get_data():
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


class BaseModelBuilder(ABC):

    @abstractmethod
    def build():
        raise NotImplementedError


class BaseModelValidator(ABC):

    @abstractmethod
    def validate():
        raise NotImplementedError
