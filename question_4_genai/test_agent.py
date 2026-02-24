"""
test_agent.py
=============
Test script for ClinicalTrialDataAgent.

Runs 3 natural-language queries against the pharmaversesdtm::ae dataset and
prints the results.  No OpenAI API key is required — the Mock LLM handles
queries when OPENAI_API_KEY is not set.

Usage
-----
    python question_4_genai/test_agent.py

    # or with a real OpenAI key:
    OPENAI_API_KEY=sk-... python question_4_genai/test_agent.py
"""

import sys
import os

# Allow running from repo root or from within question_4_genai/
sys.path.insert(0, os.path.join(os.path.dirname(__file__)))

from clinical_trial_data_agent import ClinicalTrialDataAgent

# ── Initialise agent ──────────────────────────────────────────────────────────
print("\n" + "═" * 62)
print("  ClinicalTrialDataAgent — Test Run")
print("═" * 62)

agent = ClinicalTrialDataAgent(
    csv_path=os.path.join(os.path.dirname(__file__), "adae.csv")
)
print("\n[Dataset Summary]")
print(agent.summary())

# ── Define 3 test queries ─────────────────────────────────────────────────────
TEST_QUERIES = [
    # Query 1: Severity filter (AESEV = MODERATE)
    "Give me the subjects who had Adverse Events of Moderate severity.",

    # Query 2: Body system filter (AESOC contains CARDIAC)
    "Which patients experienced Cardiac adverse events?",

    # Query 3: Specific AE term filter (AETERM contains PRURITUS)
    'Show me all subjects who reported "Pruritus" as an adverse event.',
]

# ── Run each query ────────────────────────────────────────────────────────────
print("\n" + "═" * 62)
print("  Running 3 Test Queries")
print("═" * 62 + "\n")

for i, question in enumerate(TEST_QUERIES, start=1):
    print(f"\n{'─' * 62}")
    print(f"  QUERY {i} of {len(TEST_QUERIES)}")
    print(f"{'─' * 62}")

    result = agent.ask(question)
    print(result)

    # Extra: show first few matching rows from the AE dataframe
    _, subject_ids, _ = agent.execute_query(result.intent)
    if subject_ids:
        col = result.intent.target_column
        val = result.intent.filter_value
        op  = result.intent.operation
        if op == "contains":
            mask = agent.ae[col].str.contains(val, case=False, na=False)
        else:
            mask = agent.ae[col].str.upper() == val.upper()
        sample = (
            agent.ae[mask][["USUBJID", "AETERM", "AESOC", "AESEV", "AESER", "AEOUT"]]
            .drop_duplicates()
            .head(5)
        )
        print(f"\n  Sample matching records (up to 5):\n")
        # Pretty-print with indentation
        for _, row in sample.iterrows():
            print(
                f"    USUBJID={row['USUBJID']}  AETERM={row['AETERM']!r}  "
                f"AESEV={row['AESEV']}  AESER={row['AESER']}  AEOUT={row['AEOUT']}"
            )

print("\n" + "═" * 62)
print("  All queries completed successfully.")
print("═" * 62 + "\n")
