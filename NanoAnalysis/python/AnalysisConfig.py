from pydantic import BaseModel, validator
from typing import Literal, List, Optional

class SampleSelection(BaseModel):
    years: List[int]
    origin: List[Literal['MC', 'data']]
    
    mode: Literal['name', 'process', 'all']
    by_name: List[str] = []
    by_process: List[str] = []

    @validator('years', each_item=True)
    def valid_years(cls, v):
        assert 2016 <= v <= 2025, f"Anno fuori range: {v}"
        return v

class AnalysisParameters(BaseModel):
    analyzer: str
    regions : List[Literal['4P', '3P', '2P']]
    
    
class AnalysisConfig(BaseModel):
    samples  : SampleSelection
    analysis : AnalysisParameters
