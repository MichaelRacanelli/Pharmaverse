"""
clinical_trial_data_agent.py
============================
GenAI Clinical Data Assistant — translates natural language questions into
structured Pandas queries against an Adverse Events (AE) dataset.

Architecture
------------
User Question
    │
    ▼
ClinicalTrialDataAgent.ask(question)
    │
    ▼
LangChain + OpenAI (GPT-4o-mini)
    Prompt: schema description + user question
    ↓
    Structured JSON output → QueryIntent
    │
    ▼
ClinicalTrialDataAgent.execute_query(intent)
    │
    ▼
Pandas filter on ae DataFrame
    │
    ▼
QueryResult (n_subjects, subject_ids, matched_records)

Requirements
------------
    pip install openai langchain langchain-openai pydantic pandas

Usage
-----
    from clinical_trial_data_agent import ClinicalTrialDataAgent

    agent = ClinicalTrialDataAgent("question_4_genai/adae.csv")
    result = agent.ask("Which subjects had moderate severity adverse events?")
    print(result)
"""

from __future__ import annotations

import json
import re
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd
from pydantic import BaseModel, Field
from dotenv import load_dotenv

# Load OPENAI_API_KEY (and any other vars) from .env if present
load_dotenv()

# ── Logging ──────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("ClinicalTrialDataAgent")


# ═══════════════════════════════════════════════════════════════════════════════
# 1. DATA MODELS
# ═══════════════════════════════════════════════════════════════════════════════

class QueryIntent(BaseModel):
    """
    Structured output produced by the LLM for a single user question.

    Attributes
    ----------
    target_column : str
        The AE dataset column to filter on (e.g. ``"AESEV"``).
    filter_value : str
        The value to match in ``target_column`` (e.g. ``"MODERATE"``).
    operation : str
        One of ``"equals"`` (exact match) or ``"contains"`` (substring, case-insensitive).
    reasoning : str
        Brief explanation of how the LLM mapped the question to the column/value.
    """

    target_column: str = Field(
        description="Column name in the AE dataset to filter on."
    )
    filter_value: str = Field(
        description="Value to search for in the target column."
    )
    operation: str = Field(
        default="equals",
        description='Filter operation: "equals" for exact match, "contains" for substring.'
    )
    reasoning: str = Field(
        default="",
        description="Brief explanation of the mapping decision."
    )


@dataclass
class QueryResult:
    """
    Result returned after executing a QueryIntent against the AE DataFrame.

    Attributes
    ----------
    question : str
        Original user question.
    intent : QueryIntent
        Parsed intent produced by the LLM.
    n_subjects : int
        Number of unique subjects (USUBJID) matching the filter.
    subject_ids : list[str]
        Sorted list of matching USUBJID values.
    matched_records : int
        Total number of AE records (rows) that satisfy the filter.
    """

    question: str
    intent: QueryIntent
    n_subjects: int
    subject_ids: list[str] = field(default_factory=list)
    matched_records: int = 0

    def __str__(self) -> str:  # noqa: D105
        lines = [
            "─" * 62,
            f"  Question : {self.question}",
            f"  Column   : {self.intent.target_column}",
            f"  Filter   : {self.intent.filter_value!r}  [{self.intent.operation}]",
            f"  Reasoning: {self.intent.reasoning}",
            "  ─────────────────────────────────────────────────────",
            f"  Matching AE records : {self.matched_records}",
            f"  Unique subjects (n) : {self.n_subjects}",
            f"  Subject IDs         : {self.subject_ids}",
            "─" * 62,
        ]
        return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════════
# 2. SCHEMA DEFINITION
# ═══════════════════════════════════════════════════════════════════════════════

# Human-readable data dictionary for key AE columns.
# This is injected into the LLM prompt so the model understands the schema.
AE_SCHEMA: dict[str, dict] = {
    "USUBJID": {
        "description": "Unique Subject Identifier. Composite key: study-site-subject.",
        "type": "string",
        "example": "01-701-1015",
    },
    "AETERM": {
        "description": (
            "Verbatim adverse event term as reported. "
            "Use for specific AE condition queries (e.g. 'Headache', 'Diarrhoea', "
            "'Application Site Erythema'). Case-insensitive contains match recommended."
        ),
        "type": "string",
        "example": "DIARRHOEA",
    },
    "AEDECOD": {
        "description": (
            "MedDRA preferred term (dictionary-coded AE term). "
            "More standardised than AETERM."
        ),
        "type": "string",
        "example": "DIARRHOEA",
    },
    "AESOC": {
        "description": (
            "MedDRA System Organ Class (SOC). "
            "Use when the user asks about a body system or organ (e.g. 'Cardiac', "
            "'Skin', 'Gastrointestinal', 'Nervous System'). "
            "Contains match recommended."
        ),
        "type": "string",
        "example": "CARDIAC DISORDERS",
    },
    "AEBODSYS": {
        "description": "Body System — similar to AESOC; MedDRA body system name.",
        "type": "string",
        "example": "CARDIAC DISORDERS",
    },
    "AESEV": {
        "description": (
            "AE severity / intensity. "
            "Use when user asks about severity, intensity, or grade. "
            'Valid exact values: "MILD", "MODERATE", "SEVERE".'
        ),
        "type": "categorical",
        "values": ["MILD", "MODERATE", "SEVERE"],
    },
    "AESER": {
        "description": (
            "Serious AE flag. "
            "Use when user asks about serious adverse events. "
            '"Y" = serious, "N" = not serious.'
        ),
        "type": "categorical",
        "values": ["Y", "N"],
    },
    "AEOUT": {
        "description": (
            "AE outcome. "
            "Use when user asks about outcome, resolution, or fatal events. "
            'Values: "RECOVERED/RESOLVED", "NOT RECOVERED/NOT RESOLVED", "FATAL".'
        ),
        "type": "categorical",
        "values": ["RECOVERED/RESOLVED", "NOT RECOVERED/NOT RESOLVED", "FATAL"],
    },
    "AEREL": {
        "description": (
            "Relationship to study treatment. "
            "Use when user asks about drug-related or treatment-related AEs. "
            'Values: "PROBABLE", "POSSIBLE", "REMOTE", "NONE".'
        ),
        "type": "categorical",
        "values": ["PROBABLE", "POSSIBLE", "REMOTE", "NONE"],
    },
    "AESTDTC": {
        "description": "AE start date (ISO 8601 format: YYYY-MM-DD).",
        "type": "date string",
        "example": "2012-08-05",
    },
}


def build_schema_description() -> str:
    """
    Render ``AE_SCHEMA`` as a compact, LLM-friendly text block.

    Returns
    -------
    str
        Multi-line schema description suitable for injection into a system prompt.
    """
    lines = ["ADVERSE EVENTS (AE) DATASET SCHEMA", "=" * 40]
    for col, meta in AE_SCHEMA.items():
        line = f"  {col}: {meta['description']}"
        if "values" in meta:
            line += f" | Values: {meta['values']}"
        elif "example" in meta:
            line += f" | Example: {meta['example']}"
        lines.append(line)
    return "\n".join(lines)


# ═══════════════════════════════════════════════════════════════════════════════
# 3. LANGCHAIN LLM
# ═══════════════════════════════════════════════════════════════════════════════

_SYSTEM_PROMPT = """\
You are a clinical data parsing assistant for a pharmaceutical adverse events dataset.

{schema}

Your task: Given a user question, identify the best column to filter and the value to search for.

Rules:
- target_column must be one of the column names listed in the schema above.
- filter_value must be UPPERCASE for categorical columns (AESEV, AESER, AEOUT, AEREL).
- For text search columns (AETERM, AESOC, AEDECOD), use UPPERCASE and set operation to "contains".
- For exact categorical columns (AESEV, AESER, AEOUT, AEREL), set operation to "equals".
- Provide brief, clear reasoning.

Return ONLY a valid JSON object matching this schema:
{{
  "target_column": "<column_name>",
  "filter_value": "<value_to_filter>",
  "operation": "equals" | "contains",
  "reasoning": "<brief explanation>"
}}
""".format(schema=build_schema_description())


def _langchain_parse(question: str, model: str = "gpt-4o-mini") -> QueryIntent:
    """
    Use LangChain + OpenAI to parse a question into a ``QueryIntent``.

    Parameters
    ----------
    question : str
        Natural language question from the clinical reviewer.
    model : str
        OpenAI model name (default ``"gpt-4o-mini"`` — fast and cheap).

    Returns
    -------
    QueryIntent
        Structured intent parsed from LLM JSON output.

    Raises
    ------
    ValueError
        If the LLM returns malformed JSON.
    """
    from langchain_openai import ChatOpenAI
    from langchain_core.messages import HumanMessage, SystemMessage

    llm = ChatOpenAI(model=model, temperature=0)
    messages = [
        SystemMessage(content=_SYSTEM_PROMPT),
        HumanMessage(content=question),
    ]
    response = llm.invoke(messages)
    raw = response.content.strip()

    # Strip potential markdown code fence
    raw = re.sub(r"^```(?:json)?\s*|\s*```$", "", raw, flags=re.MULTILINE).strip()

    try:
        data = json.loads(raw)
        return QueryIntent(**data)
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise ValueError(
            f"LLM returned invalid JSON:\n{raw}"
        ) from exc


# ═══════════════════════════════════════════════════════════════════════════════
# 4. MAIN AGENT CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class ClinicalTrialDataAgent:
    """
    Translates free-text clinical questions into Pandas queries on an AE dataset.

    The agent follows a three-step pipeline:
    1. **Prompt** – Inject dataset schema + user question into LLM context.
    2. **Parse**  – LLM returns structured JSON (``QueryIntent``).
    3. **Execute** – Apply the intent as a Pandas filter; return ``QueryResult``.

    Parameters
    ----------
    csv_path : str | Path
        Path to the AE CSV file (exported from ``pharmaversesdtm::ae``).
    openai_model : str
        OpenAI model to use.  Defaults to ``"gpt-4o-mini"``.
    """

    def __init__(
        self,
        csv_path: str | Path = "question_4_genai/adae.csv",
        openai_model: str = "gpt-4o-mini",
    ) -> None:
        self._csv_path = Path(csv_path)
        self._model = openai_model

        # ── Load dataset ──────────────────────────────────────────────────
        if not self._csv_path.exists():
            raise FileNotFoundError(
                f"AE dataset not found: {self._csv_path}\n"
                "Export it with: write.csv(pharmaversesdtm::ae, 'question_4_genai/adae.csv')"
            )

        self.ae: pd.DataFrame = pd.read_csv(self._csv_path, dtype=str)
        # Normalise: strip whitespace from string columns
        for col in self.ae.select_dtypes("object").columns:
            self.ae[col] = self.ae[col].str.strip()

        logger.info(
            "Loaded AE dataset: %d rows × %d cols from %s",
            len(self.ae),
            len(self.ae.columns),
            self._csv_path,
        )
        logger.info("LLM backend: OpenAI %s", self._model)

        # ── Schema ────────────────────────────────────────────────────────
        self.schema_description: str = build_schema_description()

    # ─────────────────────────────────────────────────────────────────────
    # Public API
    # ─────────────────────────────────────────────────────────────────────

    def parse_question(self, question: str) -> QueryIntent:
        """
        Send *question* to the LLM and return a ``QueryIntent``.

        Parameters
        ----------
        question : str
            Free-text natural-language question from the clinical reviewer.

        Returns
        -------
        QueryIntent
            Structured mapping: ``target_column``, ``filter_value``, ``operation``.
        """
        logger.info("Parsing question: %r", question)
        intent = _langchain_parse(question, model=self._model)
        logger.info("LLM intent: %s", intent.model_dump())
        return intent

    def execute_query(self, intent: QueryIntent) -> tuple[int, list[str], int]:
        """
        Apply *intent* as a Pandas filter on the AE DataFrame.

        Parameters
        ----------
        intent : QueryIntent
            Parsed intent from ``parse_question()``.

        Returns
        -------
        tuple[int, list[str], int]
            ``(n_subjects, subject_ids, matched_records)``
        """
        col = intent.target_column
        val = intent.filter_value

        if col not in self.ae.columns:
            logger.warning(
                "Column %r not found in dataset. Available: %s",
                col,
                list(self.ae.columns),
            )
            return 0, [], 0

        if intent.operation == "contains":
            mask = self.ae[col].str.contains(val, case=False, na=False)
        else:  # equals
            mask = self.ae[col].str.upper() == val.upper()

        filtered = self.ae[mask]
        subject_ids = sorted(filtered["USUBJID"].dropna().unique().tolist())
        return len(subject_ids), subject_ids, len(filtered)

    def ask(self, question: str) -> QueryResult:
        """
        Full pipeline: parse the question then execute the query.

        This is the primary entry point for the agent.

        Parameters
        ----------
        question : str
            Free-text question from the clinical reviewer.

        Returns
        -------
        QueryResult
            Contains the matched subject IDs, count, and reasoning.

        Examples
        --------
        >>> agent = ClinicalTrialDataAgent("question_4_genai/adae.csv")
        >>> result = agent.ask("Give me subjects who had moderate severity AEs")
        >>> print(result.n_subjects)
        """
        intent = self.parse_question(question)
        n_subjects, subject_ids, matched_records = self.execute_query(intent)
        return QueryResult(
            question=question,
            intent=intent,
            n_subjects=n_subjects,
            subject_ids=subject_ids,
            matched_records=matched_records,
        )

    # convenience ─────────────────────────────────────────────────────────
    def summary(self) -> str:
        """Return a brief text description of the loaded dataset."""
        return (
            f"Dataset  : {self._csv_path}\n"
            f"Rows     : {len(self.ae)}\n"
            f"Columns  : {list(self.ae.columns)}\n"
            f"Subjects : {self.ae['USUBJID'].nunique()}\n"
            f"LLM      : OpenAI {self._model}\n"
        )
