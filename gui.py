"""
bbsplit.gui
===========
Desktop interface (PySide6) for bbsplit.

Tabs:
  1. Input         - load CSV / TXT / SMI / SDF, or paste SMILES. Accepted
                     molecules are listed with a rendered structure preview.
  2. Rules         - enable/disable rules, load custom YAML, set depth and the
                     number of parallel workers. Changes apply immediately.
  3. Disconnections- read-only viewer: for each molecule, every disconnection,
                     showing BOTH building blocks side by side.
  4. Enumerate     - build the block pool, choose a mode + workers, run, and
                     view product structures + sortable molecular descriptors.

All chemistry lives in bbsplit.core / parallel / descriptors / io_formats.
"""

from __future__ import annotations

import os
import sys
import csv

from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFileDialog, QPlainTextEdit, QComboBox, QSpinBox,
    QListWidget, QListWidgetItem, QTableWidget, QTableWidgetItem, QCheckBox,
    QMessageBox, QSplitter, QHeaderView, QAbstractItemView, QFrame,
)
from PySide6.QtCore import Qt, QSize
from PySide6.QtGui import QPixmap, QIcon

from .core import RuleSet, SplitEngine, Disconnection, combine_two
from .parallel import split_molecules, enumerate_blocks_parallel
from .descriptors import annotate_rows, DESCRIPTOR_FIELDS
from .render import smiles_to_png
from . import io_formats


# --------------------------------------------------------------------------- #
#  Theme
# --------------------------------------------------------------------------- #
STYLESHEET = """
QWidget { font-size: 13px; color: #1f2933; }
QMainWindow, QTabWidget::pane { background: #f5f7fa; }
QTabWidget::pane { border: 1px solid #d7dee8; border-radius: 8px; top: -1px; }
QTabBar::tab {
    background: #e4e9f0; color: #3a4a5e; padding: 8px 18px; margin-right: 2px;
    border-top-left-radius: 8px; border-top-right-radius: 8px; font-weight: 600;
}
QTabBar::tab:selected { background: #2f6fed; color: white; }
QTabBar::tab:hover:!selected { background: #d2dbe8; }
QPushButton {
    background: #2f6fed; color: white; border: none; border-radius: 6px;
    padding: 7px 14px; font-weight: 600;
}
QPushButton:hover { background: #245ad1; }
QPushButton:pressed { background: #1d4cb0; }
QPushButton:disabled { background: #aab6c6; color: #eef2f7; }
QListWidget, QTableWidget, QPlainTextEdit {
    background: white; border: 1px solid #d7dee8; border-radius: 6px;
    selection-background-color: #cfe0ff; selection-color: #143a8a;
}
QHeaderView::section {
    background: #eaeef5; color: #2a3850; padding: 6px; border: none;
    border-right: 1px solid #d7dee8; font-weight: 600;
}
QComboBox, QSpinBox {
    background: white; border: 1px solid #c3cdda; border-radius: 5px; padding: 4px 6px;
}
QComboBox:focus, QSpinBox:focus { border: 1px solid #2f6fed; }
QCheckBox { spacing: 6px; }
QLabel#hint { color: #6b7787; font-style: italic; }
QLabel#section { font-size: 14px; font-weight: 700; color: #1b2a44; }
QFrame#card {
    background: white; border: 1px solid #dde4ee; border-radius: 8px;
}
"""


# --------------------------------------------------------------------------- #
#  UI helpers
# --------------------------------------------------------------------------- #
def pixmap_for(smiles: str, w: int = 240, h: int = 160) -> QPixmap:
    png = smiles_to_png(smiles, w, h)
    pm = QPixmap()
    if png:
        pm.loadFromData(png, "PNG")
    return pm


def is_trivial(frag: str) -> bool:
    heavy = sum(1 for ch in frag if ch.isalpha() and ch != "*")
    return heavy <= 1


class NumericItem(QTableWidgetItem):
    """Table item that sorts by a numeric value rather than by string."""
    def __init__(self, value):
        text = "" if value is None else (
            str(int(value)) if float(value).is_integer() else f"{value:g}")
        super().__init__(text)
        self._value = float("-inf") if value is None else float(value)
        self.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)

    def __lt__(self, other):
        if isinstance(other, NumericItem):
            return self._value < other._value
        return super().__lt__(other)


def make_enum_id(index: int, prefix: str = "ENU", width: int = 4) -> str:
    """Enumerate ID: >=3 letters + digits, e.g. ENU0001 (1-based index)."""
    return f"{prefix}{index:0{width}d}"


def max_worker_limit() -> int:
    # Allow over-subscription beyond physical cores (user's call), capped sanely.
    return max(8, (os.cpu_count() or 1) * 2)


# --------------------------------------------------------------------------- #
#  Tab 1: input (with structure previews + multi-format import)
# --------------------------------------------------------------------------- #
class InputTab(QWidget):
    THUMB = QSize(150, 110)

    def __init__(self, app_state):
        super().__init__()
        self.state = app_state
        lay = QVBoxLayout(self)

        title = QLabel("Input molecules"); title.setObjectName("section")
        lay.addWidget(title)
        info = QLabel("Load a CSV, a SMILES list (.smi/.txt), or an SD file "
                      "(.sdf/.mol) — or paste SMILES below (one per line, with an "
                      "optional name as a second token).")
        info.setObjectName("hint"); info.setWordWrap(True)
        lay.addWidget(info)

        row = QHBoxLayout()
        self.btn_csv = QPushButton("Load CSV…")
        self.btn_csv.clicked.connect(self.load_csv)
        row.addWidget(self.btn_csv)
        self.btn_file = QPushButton("Load TXT / SMI / SDF…")
        self.btn_file.clicked.connect(self.load_file)
        row.addWidget(self.btn_file)
        row.addWidget(QLabel("SMILES col:"))
        self.col_combo = QComboBox(); self.col_combo.setMinimumWidth(120)
        row.addWidget(self.col_combo)
        row.addWidget(QLabel("ID col:"))
        self.id_combo = QComboBox(); self.id_combo.setMinimumWidth(120)
        row.addWidget(self.id_combo)
        row.addStretch()
        lay.addLayout(row)

        split = QSplitter(Qt.Horizontal)

        # left: paste box + accept
        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(0, 0, 0, 0)
        ll.addWidget(QLabel("Paste SMILES:"))
        self.text = QPlainTextEdit()
        self.text.setPlaceholderText(
            "CCC(Nc1cccc(C(N)=O)n1)c1ccc(Cl)c(Cl)c1 my_compound\n…")
        ll.addWidget(self.text)
        self.btn_load = QPushButton("Accept input  →  split")
        self.btn_load.clicked.connect(self.accept_input)
        ll.addWidget(self.btn_load)
        self.status = QLabel(""); self.status.setObjectName("hint")
        ll.addWidget(self.status)
        split.addWidget(left)

        # right: accepted molecules with structure previews
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(QLabel("Accepted molecules (with structures):"))
        self.gallery = QListWidget()
        self.gallery.setViewMode(QListWidget.IconMode)
        self.gallery.setIconSize(self.THUMB)
        self.gallery.setResizeMode(QListWidget.Adjust)
        self.gallery.setMovement(QListWidget.Static)
        self.gallery.setSpacing(10)
        self.gallery.setWordWrap(True)
        rl.addWidget(self.gallery)
        split.addWidget(right)
        split.setSizes([380, 560])
        lay.addWidget(split)

        self._rows: list[dict] = []
        self._csv_path: str | None = None

    # ---- CSV (needs column choice) ----
    def load_csv(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose a CSV", "", "CSV (*.csv);;All files (*)")
        if not path:
            return
        self._csv_path = path
        cols = io_formats.csv_columns(path)
        self._rows = io_formats.read_csv  # marker; actual read on accept
        self.col_combo.clear(); self.col_combo.addItems(cols)
        self.id_combo.clear(); self.id_combo.addItems(["(none)"] + cols)
        for c in cols:
            if c.lower() in ("smiles", "smi", "structure"):
                self.col_combo.setCurrentText(c); break
        self.status.setText(f"Loaded CSV header from {os.path.basename(path)}. "
                            f"Pick the SMILES column, then Accept.")

    # ---- TXT / SMI / SDF (auto) ----
    def load_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose a file", "",
            "Molecules (*.smi *.txt *.sdf *.mol);;All files (*)")
        if not path:
            return
        try:
            items = io_formats.detect_and_read(path)
        except Exception as e:
            QMessageBox.critical(self, "Read error", str(e)); return
        if not items:
            QMessageBox.warning(self, "Empty", "No molecules read from file.")
            return
        self.state.input_molecules = items
        self._csv_path = None
        self.state.resplit_all()
        self.populate_gallery()
        self.status.setText(f"Loaded {len(items)} molecules from "
                            f"{os.path.basename(path)} and split them.")
        self.state.notify_input_ready()

    def accept_input(self):
        items: list[tuple[str, str]] = []
        if self._csv_path and self.col_combo.currentText():
            icol = self.id_combo.currentText()
            icol = None if icol in ("", "(none)") else icol
            items = io_formats.read_csv(self._csv_path,
                                        self.col_combo.currentText(), icol)
        else:
            for i, line in enumerate(self.text.toPlainText().splitlines()):
                line = line.strip()
                if not line:
                    continue
                parts = line.split(None, 1)
                smi = parts[0]
                name = parts[1].strip() if len(parts) > 1 else f"mol_{i+1}"
                items.append((name or f"mol_{i+1}", smi))
        if not items:
            QMessageBox.warning(self, "Empty", "No SMILES found.")
            return
        self.state.input_molecules = items
        self.state.resplit_all()
        self.populate_gallery()
        self.status.setText(f"Accepted {len(items)} molecules and split them. "
                            f"See the 'Disconnections' tab.")
        self.state.notify_input_ready()

    def populate_gallery(self):
        self.gallery.clear()
        from rdkit import Chem
        for mid, smi in self.state.input_molecules:
            valid = Chem.MolFromSmiles(smi) is not None
            it = QListWidgetItem(f"{mid}\n{smi if len(smi) < 40 else smi[:37]+'…'}")
            it.setTextAlignment(Qt.AlignHCenter | Qt.AlignBottom)
            if valid:
                it.setIcon(QIcon(pixmap_for(smi, self.THUMB.width(),
                                            self.THUMB.height())))
            else:
                it.setText(f"{mid}\n(invalid SMILES)")
            it.setSizeHint(QSize(self.THUMB.width() + 30, self.THUMB.height() + 44))
            self.gallery.addItem(it)


# --------------------------------------------------------------------------- #
#  Tab 2: rules (changes apply immediately)
# --------------------------------------------------------------------------- #
class RulesTab(QWidget):
    def __init__(self, app_state):
        super().__init__()
        self.state = app_state
        self._building = False
        lay = QVBoxLayout(self)

        title = QLabel("Disconnection rules"); title.setObjectName("section")
        lay.addWidget(title)
        note = QLabel("Any change here is applied immediately: all input "
                      "molecules are re-split with the currently enabled rules.")
        note.setObjectName("hint"); note.setWordWrap(True)
        lay.addWidget(note)

        top = QHBoxLayout()
        top.addWidget(QLabel("Depth (max bonds):"))
        self.depth = QSpinBox(); self.depth.setRange(1, 3); self.depth.setValue(2)
        self.depth.valueChanged.connect(self.on_depth)
        top.addWidget(self.depth)
        top.addWidget(QLabel("Parallel workers:"))
        self.workers = QSpinBox()
        self.workers.setRange(1, max_worker_limit())
        self.workers.setValue(min(4, max_worker_limit()))
        self.workers.setToolTip("Worker processes for splitting. May exceed core "
                                "count if you wish (each runs as its own process).")
        self.workers.valueChanged.connect(self.on_workers)
        top.addWidget(self.workers)
        self.btn_yaml = QPushButton("Add rules from YAML…")
        self.btn_yaml.clicked.connect(self.load_yaml)
        top.addWidget(self.btn_yaml)
        top.addStretch()
        lay.addLayout(top)

        self.list = QListWidget()
        self.list.itemChanged.connect(self.on_toggle)
        lay.addWidget(self.list)

        self.status = QLabel(""); self.status.setObjectName("hint")
        lay.addWidget(self.status)
        self.refresh_list()

    def refresh_list(self):
        self._building = True
        self.list.clear()
        for r in self.state.ruleset.rules:
            it = QListWidgetItem(f"[{r.priority}] {r.name} — {r.description}")
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            it.setCheckState(Qt.Checked if r.enabled else Qt.Unchecked)
            it.setData(Qt.UserRole, r.name)
            self.list.addItem(it)
        self._building = False

    def on_toggle(self, item):
        if self._building:
            return
        name = item.data(Qt.UserRole)
        for r in self.state.ruleset.rules:
            if r.name == name:
                r.enabled = (item.checkState() == Qt.Checked)
        self._apply()

    def on_depth(self, v):
        self.state.engine.max_bonds = v
        self._apply()

    def on_workers(self, v):
        self.state.workers = v
        # No re-split needed just for worker count, but keep it cheap & correct.
        self._apply()

    def _apply(self):
        n = self.state.resplit_all()
        if n is not None:
            self.status.setText(f"Applied. Re-split {n} molecules "
                                f"(workers={self.state.workers}).")
        self.state.notify_rules_applied()

    def load_yaml(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Rules YAML", "", "YAML (*.yaml *.yml)")
        if not path:
            return
        try:
            self.state.ruleset.merge_yaml(path, override=True)
        except Exception as e:
            QMessageBox.critical(self, "YAML error", str(e)); return
        self.refresh_list()
        self._apply()
        QMessageBox.information(self, "OK", "Rules added/updated and applied.")


# --------------------------------------------------------------------------- #
#  Tab 3: disconnections (read-only; shows BOTH blocks)
# --------------------------------------------------------------------------- #
class SplitTab(QWidget):
    BLK = QSize(260, 175)

    def __init__(self, app_state):
        super().__init__()
        self.state = app_state
        outer = QVBoxLayout(self)
        title = QLabel("Disconnections"); title.setObjectName("section")
        outer.addWidget(title)

        lay = QHBoxLayout()
        # left - molecule list
        left = QVBoxLayout()
        left.addWidget(QLabel("Molecules:"))
        self.mol_list = QListWidget()
        self.mol_list.currentRowChanged.connect(self.show_disconnections)
        left.addWidget(self.mol_list)
        lw = QWidget(); lw.setLayout(left)

        # middle - list of ways
        mid = QVBoxLayout()
        mid.addWidget(QLabel("Ways to disconnect (read-only):"))
        self.disc_list = QListWidget()
        self.disc_list.currentItemChanged.connect(self._preview)
        mid.addWidget(self.disc_list)
        mw = QWidget(); mw.setLayout(mid)

        # right - the resulting building blocks, side by side
        right = QVBoxLayout()
        right.addWidget(QLabel("Resulting building blocks:"))
        self.blocks_row = QHBoxLayout()
        self.block_frames: list[tuple[QLabel, QLabel]] = []
        blocks_holder = QWidget(); blocks_holder.setLayout(self.blocks_row)
        right.addWidget(blocks_holder)
        right.addStretch()
        rw = QWidget(); rw.setLayout(right)

        sp = QSplitter(Qt.Horizontal)
        sp.addWidget(lw); sp.addWidget(mw); sp.addWidget(rw)
        sp.setSizes([240, 360, 560])
        lay.addWidget(sp)
        outer.addLayout(lay)

    def refresh_molecules(self):
        self.mol_list.clear()
        for mid, smi in self.state.input_molecules:
            n = len(self.state.disconnections.get(mid, []))
            self.mol_list.addItem(f"{mid}  ({n} ways)")
        if self.mol_list.count():
            self.mol_list.setCurrentRow(0)

    def show_disconnections(self, row):
        self.disc_list.clear()
        self._clear_blocks()
        if row < 0 or row >= len(self.state.input_molecules):
            return
        mid, _ = self.state.input_molecules[row]
        for i, d in enumerate(self.state.disconnections.get(mid, [])):
            it = QListWidgetItem(f"#{i}  ·  {d.n_blocks} blocks  ·  {d.label()}")
            it.setData(Qt.UserRole, (mid, i))
            self.disc_list.addItem(it)
        if self.disc_list.count():
            self.disc_list.setCurrentRow(0)

    def _clear_blocks(self):
        while self.blocks_row.count():
            item = self.blocks_row.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()
        self.block_frames = []

    def _preview(self, cur, _prev):
        self._clear_blocks()
        if not cur:
            return
        mid, i = cur.data(Qt.UserRole)
        d = self.state.disconnections[mid][i]
        # show EVERY resulting block side by side (2 or 3 blocks)
        for k, frag in enumerate(d.fragments):
            card = QFrame(); card.setObjectName("card")
            v = QVBoxLayout(card)
            cap = QLabel(f"Block {chr(ord('A')+k)}")
            cap.setAlignment(Qt.AlignHCenter)
            cap.setStyleSheet("font-weight:700; color:#2f6fed;")
            v.addWidget(cap)
            img = QLabel(); img.setAlignment(Qt.AlignCenter)
            img.setMinimumSize(self.BLK)
            pm = pixmap_for(frag, self.BLK.width(), self.BLK.height())
            if pm.isNull():
                img.setText("(cannot render)")
            else:
                img.setPixmap(pm)
            v.addWidget(img)
            smi = QLabel(frag); smi.setWordWrap(True)
            smi.setAlignment(Qt.AlignHCenter); smi.setObjectName("hint")
            v.addWidget(smi)
            self.blocks_row.addWidget(card)
        self.blocks_row.addStretch()


# --------------------------------------------------------------------------- #
#  Tab 4: enumerate
# --------------------------------------------------------------------------- #
class EnumTab(QWidget):
    THUMB = QSize(220, 150)

    def __init__(self, app_state):
        super().__init__()
        self.state = app_state
        lay = QVBoxLayout(self)
        title = QLabel("Enumerate products"); title.setObjectName("section")
        lay.addWidget(title)

        opt = QHBoxLayout()
        self.dedup = QCheckBox("Deduplicate products"); self.dedup.setChecked(True)
        opt.addWidget(self.dedup)
        self.skip_trivial = QCheckBox("Skip trivial blocks ([*]N, etc.)")
        self.skip_trivial.setChecked(True)
        opt.addWidget(self.skip_trivial)
        opt.addStretch()
        self.btn_pool = QPushButton("Build block pool")
        self.btn_pool.clicked.connect(self.build_pool)
        opt.addWidget(self.btn_pool)
        lay.addLayout(opt)

        self.pool_info = QLabel("Block pool not built yet.")
        self.pool_info.setObjectName("hint")
        lay.addWidget(self.pool_info)

        mode = QHBoxLayout()
        mode.addWidget(QLabel("Mode:"))
        self.mode = QComboBox()
        self.mode.addItems([
            "All × all (full Cartesian, 2 roles)",
            "Reproduce input molecules only",
            "One chosen block A × all B",
        ])
        mode.addWidget(self.mode)
        mode.addWidget(QLabel("Fixed block:"))
        self.role_combo = QComboBox(); self.role_combo.setMinimumWidth(160)
        mode.addWidget(self.role_combo)
        mode.addWidget(QLabel("Workers:"))
        self.workers = QSpinBox()
        self.workers.setRange(1, max_worker_limit())
        self.workers.setValue(min(4, max_worker_limit()))
        self.workers.setToolTip("Worker processes for enumeration.")
        mode.addWidget(self.workers)
        mode.addStretch()
        self.btn_enum = QPushButton("Enumerate")
        self.btn_enum.clicked.connect(self.run_enum)
        mode.addWidget(self.btn_enum)
        lay.addLayout(mode)

        self.base_cols = ["ID", "Structure", "product SMILES"]
        self.all_cols = self.base_cols + DESCRIPTOR_FIELDS
        self.table = QTableWidget(0, len(self.all_cols))
        self.table.setHorizontalHeaderLabels(self.all_cols)
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.Stretch)
        self.table.setColumnWidth(0, 72)
        self.table.setColumnWidth(1, self.THUMB.width() + 12)
        self.table.setIconSize(self.THUMB)
        self.table.verticalHeader().setDefaultSectionSize(self.THUMB.height() + 8)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setAlternatingRowColors(True)
        self.table.setSortingEnabled(True)
        lay.addWidget(self.table)

        hint = QLabel("Tip: click any descriptor column header to sort the "
                      "products by that property.")
        hint.setObjectName("hint")
        lay.addWidget(hint)

        self.btn_export = QPushButton("Export CSV…")
        self.btn_export.clicked.connect(self.export)
        lay.addWidget(self.btn_export)

        self.pool_a: list[str] = []
        self.pool_b: list[str] = []
        self.results: list[dict] = []

    def build_pool(self):
        from rdkit import Chem
        a_set, b_set = set(), set()
        for mid, _ in self.state.input_molecules:
            for d in self.state.disconnections.get(mid, []):
                if d.n_blocks != 2:
                    continue
                for frag in d.fragments:
                    if self.skip_trivial.isChecked() and is_trivial(frag):
                        continue
                    m = Chem.MolFromSmiles(frag)
                    if m is None:
                        continue
                    el = None
                    for at in m.GetAtoms():
                        if at.GetAtomicNum() == 0:
                            nb = at.GetNeighbors()
                            if nb:
                                el = nb[0].GetSymbol()
                            break
                    (a_set if el in ("N", "O") else b_set).add(frag)
        self.pool_a = sorted(a_set)
        self.pool_b = sorted(b_set)
        self.role_combo.clear()
        self.role_combo.addItems(self.pool_a)
        self.pool_info.setText(
            f"Pool built: A blocks (amine/ether) = {len(self.pool_a)}, "
            f"B blocks (carbon) = {len(self.pool_b)}.")

    def run_enum(self):
        if not self.pool_a or not self.pool_b:
            QMessageBox.warning(self, "No pool", "Build the block pool first.")
            return
        dedup = self.dedup.isChecked()
        workers = self.workers.value()
        mode = self.mode.currentIndex()
        if mode == 0:
            rows = enumerate_blocks_parallel(
                self.pool_a, self.pool_b, dedup=dedup, workers=workers)
        elif mode == 2:
            fixed = self.role_combo.currentText()
            rows = enumerate_blocks_parallel(
                [fixed], self.pool_b, dedup=dedup, workers=workers)
        else:
            rows, seen = [], set()
            for mid, _ in self.state.input_molecules:
                for d in self.state.disconnections.get(mid, []):
                    if d.n_blocks != 2:
                        continue
                    p = combine_two(d.fragments[0], d.fragments[1])
                    if p and (not dedup or p not in seen):
                        seen.add(p)
                        rows.append({"product": p})
        for i, r in enumerate(rows, start=1):
            r["ID"] = make_enum_id(i)
        annotate_rows(rows, smiles_key="product")
        self.results = rows
        self._fill_table(rows)
        QMessageBox.information(self, "Done", f"Enumerated {len(rows)} products.")

    def _fill_table(self, rows):
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(rows))
        for i, r in enumerate(rows):
            self.table.setItem(i, 0, QTableWidgetItem(r.get("ID", "")))
            pm = pixmap_for(r["product"], self.THUMB.width(), self.THUMB.height())
            icon_item = QTableWidgetItem()
            if not pm.isNull():
                icon_item.setIcon(QIcon(pm))
            self.table.setItem(i, 1, icon_item)
            self.table.setItem(i, 2, QTableWidgetItem(r["product"]))
            for j, field in enumerate(DESCRIPTOR_FIELDS, start=len(self.base_cols)):
                self.table.setItem(i, j, NumericItem(r.get(field)))
        self.table.setSortingEnabled(True)

    def export(self):
        if not self.results:
            QMessageBox.warning(self, "Empty", "Nothing to export.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save CSV", "enumerated.csv", "CSV (*.csv)")
        if not path:
            return
        preferred = (["ID", "blockA", "blockB", "a_idx", "b_idx", "product"]
                     + DESCRIPTOR_FIELDS)
        present = {k for r in self.results for k in r.keys()}
        keys = [k for k in preferred if k in present]
        keys += [k for k in sorted(present) if k not in keys]
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            w.writerows(self.results)
        QMessageBox.information(self, "Saved", f"Wrote {len(self.results)} rows "
                                               f"to {path}.")


# --------------------------------------------------------------------------- #
#  Application state / main window
# --------------------------------------------------------------------------- #
class AppState:
    def __init__(self):
        self.ruleset = RuleSet.default()
        self.engine = SplitEngine(self.ruleset, max_bonds=2)
        self.input_molecules: list[tuple[str, str]] = []
        self.disconnections: dict[str, list[Disconnection]] = {}
        self.workers: int = min(4, max_worker_limit())
        self._split_tab = None
        self._enum_tab = None

    def resplit_all(self):
        if not self.input_molecules:
            self.disconnections = {}
            return None
        self.disconnections = split_molecules(
            self.input_molecules, self.ruleset,
            max_bonds=self.engine.max_bonds, workers=self.workers)
        return len(self.disconnections)

    def notify_input_ready(self):
        if self._split_tab:
            self._split_tab.refresh_molecules()

    def notify_rules_applied(self):
        if self._split_tab:
            self._split_tab.refresh_molecules()


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("bbsplit — Building Block Splitter & Enumerator")
        self.resize(1120, 760)
        self.state = AppState()

        tabs = QTabWidget()
        self.input_tab = InputTab(self.state)
        self.rules_tab = RulesTab(self.state)
        self.split_tab = SplitTab(self.state)
        self.enum_tab = EnumTab(self.state)
        self.state._split_tab = self.split_tab
        self.state._enum_tab = self.enum_tab

        tabs.addTab(self.input_tab, "1. Input")
        tabs.addTab(self.rules_tab, "2. Rules")
        tabs.addTab(self.split_tab, "3. Disconnections")
        tabs.addTab(self.enum_tab, "4. Enumerate")
        self.setCentralWidget(tabs)


def main():
    app = QApplication(sys.argv)
    app.setStyleSheet(STYLESHEET)
    win = MainWindow()
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
