import os
import tempfile
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

# Use SQLite by default for local use, but keep mutable runtime data outside
# the installed application directory. Corporate application controls can make
# that directory read-only, which otherwise turns a valid assessment into an
# HTTP 500 when its audit record is persisted. DATABASE_URL still takes
# precedence for managed Postgres or an approved persistent SQLite location.
_runtime_database = Path(os.environ.get("PRISM_RUNTIME_DB") or tempfile.gettempdir()) / "PRISM" / "runtime" / "chem_engine.db"
_runtime_database.parent.mkdir(parents=True, exist_ok=True)
SQLALCHEMY_DATABASE_URL = os.getenv("DATABASE_URL", f"sqlite:///{_runtime_database.as_posix()}")

# SQLite needs check_same_thread=False, Postgres does not
connect_args = {"check_same_thread": False} if SQLALCHEMY_DATABASE_URL.startswith("sqlite") else {}

engine = create_engine(
    SQLALCHEMY_DATABASE_URL, connect_args=connect_args
)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)

Base = declarative_base()

def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()
