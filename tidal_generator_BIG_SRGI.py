"""
=============================================================
SYNTHETIC TIDAL TIME SERIES GENERATOR
Berbasis Data Aktual BIG SRGI (Prediksi Pasut)
Metode: Analisis Harmonik (utide) + Rekonstruksi

Referensi:
  Codiga, D.L. (2011). Unified Tidal Analysis and Prediction
  Using the UTide Matlab Functions. Technical Report 2011-01.
  Graduate School of Oceanography, University of Rhode Island.
  doi:10.13140/RG.2.1.3761.2008

Dikembangkan untuk penelitian pemodelan banjir pesisir
Maluku Utara — UMMU Ternate
=============================================================
"""

import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import re
import os
import csv
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import matplotlib.gridspec as gridspec

try:
    import utide
    UTIDE_OK = True
except ImportError:
    UTIDE_OK = False


# ─────────────────────────────────────────────
# WARNA & STYLE
# ─────────────────────────────────────────────
BG        = "#F4F6FA"   # latar jendela
BG2       = "#FFFFFF"   # kartu, area grafik
BG3       = "#E8EDF5"   # kolom isian, tab, bar status
ACCENT    = "#1F6FD1"   # biru utama
ACCENT2   = "#2E7D32"   # hijau (status berhasil)
WARN      = "#B36B00"   # oranye tua (hasil/peringatan)
ERR       = "#C62828"
TXT       = "#1F2937"   # teks utama
TXT2      = "#5B6474"   # teks sekunder
BORDER    = "#CBD3E1"
BTN_NEUTRAL = "#5B6B82" # tombol netral
HAT_CLR   = "#D32F2F"
LAT_CLR   = "#D32F2F"
MSL_CLR   = "#8A94A6"
TIDE_CLR  = "#1F6FD1"


def egm2008_undulation(lat, lon):
    """Undulasi geoid EGM2008 (m) di titik (lat, lon). Butuh internet.
    1) pyproj + grid PROJ (diunduh otomatis dari cdn.proj.org)
    2) cadangan: kalkulator GeoidEval GeographicLib
    Mengembalikan (N, sumber); melempar RuntimeError bila keduanya gagal."""
    errors = []
    try:
        from pyproj import Transformer, network
        network.set_network_enabled(True)
        t = Transformer.from_crs("EPSG:4979", "EPSG:4326+3855", always_xy=True)
        _, _, H = t.transform(lon, lat, 0.0)
        if np.isfinite(H) and H != 0.0 and abs(H) < 200:
            return float(-H), "pyproj/PROJ (grid EGM2008)"
        errors.append("pyproj: grid EGM2008 tidak tersedia")
    except Exception as e:
        errors.append(f"pyproj: {e}")
    try:
        import urllib.request
        url = ("https://geographiclib.sourceforge.io/cgi-bin/GeoidEval"
               f"?input={lat:.6f}+{lon:.6f}&option=Submit")
        with urllib.request.urlopen(url, timeout=20) as r:
            html = r.read().decode("utf-8", errors="ignore")
        text = re.sub(r"<[^>]+>", " ", html)
        m = re.search(r"EGM2008\s*=\s*([-\d.]+)", text)
        if m:
            return float(m.group(1)), "GeographicLib GeoidEval"
        errors.append("GeographicLib: format halaman tidak dikenali")
    except Exception as e:
        errors.append(f"GeographicLib: {e}")
    raise RuntimeError("; ".join(errors))


def compute_metrics(obs, mod):
    """Statistik validasi: obs = data BIG (m, MSL), mod = rekonstruksi (m, MSL)."""
    obs = np.asarray(obs, dtype=float); mod = np.asarray(mod, dtype=float)
    err = mod - obs
    rmse = float(np.sqrt(np.mean(err**2)))
    mae  = float(np.mean(np.abs(err)))
    bias = float(np.mean(err))
    r    = float(np.corrcoef(obs, mod)[0, 1]) if len(obs) > 2 else float('nan')
    nse  = float(1 - np.sum(err**2) / np.sum((obs - obs.mean())**2))
    return {'n': int(len(obs)), 'RMSE': rmse, 'MAE': mae, 'Bias': bias,
            'R': r, 'NSE': nse, 'max_err': float(np.max(np.abs(err)))}


def compute_validation(times_utc, z, coef, z0, start_utc, end_utc):
    """Bandingkan rekonstruksi utide dengan data BIG pada jam yang sama.
    Mengembalikan statistik setahun penuh dan statistik jendela simulasi."""
    t_np = np.array([np.datetime64(t) for t in times_utc])
    obs  = np.asarray(z, dtype=float) - z0
    mod  = np.asarray(utide.reconstruct(t_np, coef, verbose=False).h, dtype=float) - z0
    full = compute_metrics(obs, mod)
    mask = np.array([start_utc <= t <= end_utc for t in times_utc])
    win  = compute_metrics(obs[mask], mod[mask]) if mask.sum() >= 3 else None
    if win is not None:
        ow, mw = obs[mask], mod[mask]
        win['dHigh'] = float(mw.max() - ow.max())   # selisih pasang tertinggi
        win['dLow']  = float(mw.min() - ow.min())   # selisih surut terendah
    return {'full': full, 'win': win, 'obs': obs, 'mod': mod, 'mask': mask}


class TidalApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Synthetic Tidal Time Series Generator — BIG SRGI")
        self.geometry("1200x820")
        self.minsize(1000, 700)
        self.configure(bg=BG)

        # State
        self.raw_data       = []        # [(datetime, z), ...]
        self.coef           = None      # utide coefficients
        self.constituents   = {}        # {name: {A, g, T}}
        self.scaled_values  = []        # final time series
        self.time_labels    = []
        self.hat_val        = tk.DoubleVar(value=1.019)
        self.lat_val        = tk.DoubleVar(value=-0.943)
        self.msl_val        = tk.DoubleVar(value=0.000)
        self.lat_site       = tk.DoubleVar(value=2.0160)
        self.utc_offset     = tk.DoubleVar(value=9.0)    # WIT = UTC+9
        self.dem_offset     = tk.DoubleVar(value=0.0)    # koreksi MSL -> datum DEMNAS (m)
        self.lon_site       = tk.DoubleVar(value=0.0)
        self.msl_portal     = tk.DoubleVar(value=0.0)    # MSL dari Model Pasut SRGI (thd geoid)
        self.n_inageoid     = tk.DoubleVar(value=0.0)    # N dari Model Geoid SRGI
        self.n_egm2008      = tk.StringVar(value="")     # diisi otomatis
        self.z0             = 0.0                        # muka air rata-rata (datum data BIG)
        self.val            = None                       # hasil validasi
        self.start_utc      = None
        self.end_utc        = None
        self.duration_h     = tk.IntVar(value=30)
        self.interval_min   = tk.IntVar(value=15)
        self.start_dt_str   = tk.StringVar(value="2025-01-01 00:00")
        self.status_var     = tk.StringVar(value="Siap. Mulai dengan memuat file BIG SRGI.")
        self.progress_var   = tk.DoubleVar(value=0)

        self._build_ui()

    # ══════════════════════════════════════════
    # UI BUILDER
    # ══════════════════════════════════════════
    def _build_ui(self):
        # ── Top bar ──
        topbar = tk.Frame(self, bg=BG, pady=8)
        topbar.pack(fill=tk.X, padx=16)
        tk.Label(topbar, text="🌊  Synthetic Tidal Time Series Generator",
                 font=("Segoe UI", 15, "bold"), bg=BG, fg=TXT).pack(side=tk.LEFT)
        tk.Label(topbar, text="Berbasis Data Prediksi Pasut BIG SRGI  |  Metode: Analisis Harmonik (utide)",
                 font=("Segoe UI", 9), bg=BG, fg=TXT2).pack(side=tk.LEFT, padx=16)

        # ── Main paned ──
        main = tk.PanedWindow(self, orient=tk.HORIZONTAL, bg=BG,
                              sashwidth=6, sashrelief=tk.FLAT)
        main.pack(fill=tk.BOTH, expand=True, padx=10, pady=(0,8))

        left  = self._build_left(main)
        right = self._build_right(main)
        main.add(left,  minsize=320, width=370)
        main.add(right, minsize=500)

        # ── Status bar ──
        sbar = tk.Frame(self, bg=BG3, pady=4)
        sbar.pack(fill=tk.X, side=tk.BOTTOM)
        self.prog_bar = ttk.Progressbar(sbar, variable=self.progress_var,
                                        maximum=100, length=160, mode='determinate')
        self.prog_bar.pack(side=tk.LEFT, padx=10)
        tk.Label(sbar, textvariable=self.status_var, font=("Segoe UI", 9),
                 bg=BG3, fg=TXT2, anchor=tk.W).pack(side=tk.LEFT, fill=tk.X)

        # Style ttk
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TProgressbar', troughcolor=BG, background=ACCENT,
                        bordercolor=BG3, lightcolor=ACCENT, darkcolor=ACCENT)
        style.configure('TNotebook', background=BG, borderwidth=0)
        style.configure('TNotebook.Tab', background=BG3, foreground=TXT2,
                        padding=[14,6], font=('Segoe UI', 9))
        style.map('TNotebook.Tab', background=[('selected', BG2)],
                  foreground=[('selected', TXT)])

    def _card(self, parent, title, **kw):
        outer = tk.Frame(parent, bg=BG2, bd=0, relief=tk.FLAT,
                         highlightbackground=BORDER, highlightthickness=1)
        outer.pack(fill=tk.X, pady=(0,10), **kw)
        tk.Label(outer, text=title, font=("Segoe UI", 9, "bold"),
                 bg=BG2, fg=ACCENT).pack(anchor=tk.W, padx=12, pady=(8,2))
        inner = tk.Frame(outer, bg=BG2)
        inner.pack(fill=tk.X, padx=12, pady=(0,10))
        return inner

    def _row(self, parent, label, widget_fn, **kw):
        row = tk.Frame(parent, bg=BG2)
        row.pack(fill=tk.X, pady=2)
        tk.Label(row, text=label, font=("Segoe UI", 9), bg=BG2,
                 fg=TXT2, width=20, anchor=tk.W).pack(side=tk.LEFT)
        w = widget_fn(row, **kw)
        w.pack(side=tk.LEFT, fill=tk.X, expand=True)
        return w

    def _entry(self, parent, textvariable, **kw):
        e = tk.Entry(parent, textvariable=textvariable,
                     bg=BG3, fg=TXT, insertbackground=TXT,
                     relief=tk.FLAT, font=("Segoe UI", 9),
                     highlightbackground=BORDER, highlightthickness=1, **kw)
        return e

    def _btn(self, parent, text, cmd, color=ACCENT, **kw):
        b = tk.Button(parent, text=text, command=cmd,
                      bg=color, fg="white", relief=tk.FLAT,
                      font=("Segoe UI", 9, "bold"),
                      activebackground=color, activeforeground="white",
                      cursor="hand2", padx=10, pady=4, **kw)
        return b

    # ── LEFT PANEL ──
    def _build_left(self, parent):
        frame = tk.Frame(parent, bg=BG)
        canvas = tk.Canvas(frame, bg=BG, highlightthickness=0)
        scroll = tk.Scrollbar(frame, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=scroll.set)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        inner = tk.Frame(canvas, bg=BG)
        win = canvas.create_window((0,0), window=inner, anchor='nw')
        inner.bind('<Configure>', lambda e: canvas.configure(
            scrollregion=canvas.bbox('all')))
        canvas.bind('<Configure>', lambda e: canvas.itemconfig(win, width=e.width))

        # ── 1. Load file ──
        c1 = self._card(inner, "📂  1. Muat File BIG SRGI", padx=8)
        tk.Label(c1, text="File .txt atau .csv dari menu\n'Unduh Prediksi Pasut' SRGI BIG",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2, justify=tk.LEFT).pack(anchor=tk.W, pady=(0,6))
        self._btn(c1, "Pilih File BIG SRGI...", self._load_file).pack(fill=tk.X, pady=2)
        self.file_label = tk.Label(c1, text="Belum ada file dipilih",
                                   font=("Segoe UI", 8), bg=BG2, fg=TXT2,
                                   wraplength=280, justify=tk.LEFT)
        self.file_label.pack(anchor=tk.W, pady=2)
        self.data_info = tk.Label(c1, text="", font=("Segoe UI", 8, "italic"),
                                  bg=BG2, fg=ACCENT2)
        self.data_info.pack(anchor=tk.W)

        # ── 2. Parameter lokasi ──
        c2 = self._card(inner, "📍  2. Parameter Lokasi & Datum", padx=8)
        self._row(c2, "Latitude (°)", self._entry, textvariable=self.lat_site)
        self._row(c2, "Longitude (°)", self._entry, textvariable=self.lon_site)
        self._row(c2, "HAT (m, MSL ref)", self._entry, textvariable=self.hat_val)
        self._row(c2, "LAT (m, MSL ref)", self._entry, textvariable=self.lat_val)
        self._row(c2, "Zona waktu (UTC+jam)", self._entry, textvariable=self.utc_offset)
        tk.Label(c2, text="Lat/lon terisi otomatis dari file BIG. HAT, LAT dari\n"
                          "Model Pasut SRGI. WIT (Ternate) = 9, WITA = 8, WIB = 7",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2, justify=tk.LEFT).pack(anchor=tk.W, pady=(4,6))

        tk.Label(c2, text="Koreksi datum ke DEMNAS (EGM2008)",
                 font=("Segoe UI", 8, "bold"), bg=BG2, fg=ACCENT).pack(anchor=tk.W)
        self._row(c2, "MSL portal SRGI (m)", self._entry, textvariable=self.msl_portal)
        self._row(c2, "N InaGeoid2020 (m)", self._entry, textvariable=self.n_inageoid)
        self._row(c2, "N EGM2008 (m)", self._entry, textvariable=self.n_egm2008)
        self._btn(c2, "Hitung Koreksi ke DEMNAS", self._calc_dem,
                  color="#7B5EA7").pack(fill=tk.X, pady=(4,2))
        self._row(c2, "Koreksi ke DEM (m)", self._entry, textvariable=self.dem_offset)
        self.dem_info = tk.Label(c2, text="",
                                 font=("Segoe UI", 8), bg=BG2, fg=TXT2,
                                 justify=tk.LEFT, wraplength=300)
        self.dem_info.pack(anchor=tk.W, pady=(2,0))

        # ── 3. Analisis harmonik ──
        c3 = self._card(inner, "🔬  3. Analisis Harmonik (utide)", padx=8)
        self._btn(c3, "▶  Jalankan Analisis Harmonik",
                  self._run_analysis, color=ACCENT2).pack(fill=tk.X, pady=2)
        self.const_frame = tk.Frame(c3, bg=BG2)
        self.const_frame.pack(fill=tk.X, pady=4)

        # ── 4. Analisis Waktu Simulasi ──
        c4 = self._card(inner, "🔍  4. Analisis Waktu Simulasi", padx=8)
        tk.Label(c4, text="Cari tanggal terbaik untuk simulasi\nberbasis data BIG SRGI yang dimuat",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2, justify=tk.LEFT).pack(anchor=tk.W, pady=(0,4))

        # Input rising limb BC hulu
        row_rl = tk.Frame(c4, bg=BG2)
        row_rl.pack(fill=tk.X, pady=(0,6))
        tk.Label(row_rl, text="Rising limb BC hulu (jam):",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2, anchor=tk.W).pack(side=tk.LEFT)
        self.rising_limb_h = tk.DoubleVar(value=11.0)
        e_rl = tk.Entry(row_rl, textvariable=self.rising_limb_h, width=7,
                        bg=BG3, fg=TXT, insertbackground=TXT, relief=tk.FLAT,
                        font=("Segoe UI", 9),
                        highlightbackground="#7B5EA7", highlightthickness=1)
        e_rl.pack(side=tk.LEFT, padx=(6,0))
        tk.Label(row_rl, text="jam  (dari HEC-HMS)",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2).pack(side=tk.LEFT, padx=4)

        self._btn(c4, "🌕  Cari Spring Tide (z tertinggi)",
                  self._find_spring_tide, color="#7B5EA7").pack(fill=tk.X, pady=2)
        self._btn(c4, "📅  Tampilkan Top 10 z Tertinggi",
                  self._show_top_z, color="#5A7A8A").pack(fill=tk.X, pady=2)
        self.spring_result = tk.Label(c4, text="", font=("Segoe UI", 8),
                                      bg=BG2, fg=WARN, wraplength=280, justify=tk.LEFT)
        self.spring_result.pack(anchor=tk.W, pady=3)

        # ── 5. Parameter Generate ──
        c5 = self._card(inner, "⚙️  5. Parameter Generate", padx=8)
        self._row(c5, "Waktu mulai", self._entry, textvariable=self.start_dt_str)
        self._row(c5, "Durasi (jam)", self._entry, textvariable=self.duration_h)
        self._row(c5, "Interval (menit)", self._entry, textvariable=self.interval_min)
        tk.Label(c5, text="Format waktu: YYYY-MM-DD HH:MM",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2).pack(anchor=tk.W, pady=(2,0))
        self._btn(c5, "🌊  Generate Time Series",
                  self._generate, color=ACCENT).pack(fill=tk.X, pady=(4,2))

        # ── 6. Validasi ──
        cv = self._card(inner, "✅  6. Validasi terhadap Data BIG", padx=8)
        tk.Label(cv, text="Bandingkan rekonstruksi dengan data per jam BIG\n"
                          "(jendela simulasi dan setahun penuh)",
                 font=("Segoe UI", 8), bg=BG2, fg=TXT2, justify=tk.LEFT).pack(anchor=tk.W, pady=(0,4))
        self._btn(cv, "▶  Jalankan Validasi", self._run_validation,
                  color="#3E8E7E").pack(fill=tk.X, pady=2)
        self.val_info = tk.Label(cv, text="", font=("Consolas", 8), bg=BG2,
                                 fg=TXT, justify=tk.LEFT, anchor=tk.W)
        self.val_info.pack(anchor=tk.W, pady=3)

        # ── 7. Export ──
        c6 = self._card(inner, "💾  7. Export", padx=8)
        self._btn(c6, "Export CSV", self._export_csv).pack(fill=tk.X, pady=2)
        self._btn(c6, "Export Grafik (PNG)", self._export_png,
                  color=BTN_NEUTRAL).pack(fill=tk.X, pady=2)
        self._btn(c6, "Export Validasi (CSV + PNG)", self._export_validation,
                  color="#3E8E7E").pack(fill=tk.X, pady=2)
        self.export_info = tk.Label(c6, text="", font=("Segoe UI", 8),
                                    bg=BG2, fg=ACCENT2, wraplength=280)
        self.export_info.pack(anchor=tk.W)

        return frame

    # ── RIGHT PANEL ──
    def _build_right(self, parent):
        frame = tk.Frame(parent, bg=BG)
        nb = ttk.Notebook(frame)
        nb.pack(fill=tk.BOTH, expand=True, padx=8, pady=8)

        # Tab 1 — Grafik
        tab_plot = tk.Frame(nb, bg=BG)
        nb.add(tab_plot, text="  📈  Grafik Time Series  ")
        self._build_plot_tab(tab_plot)

        # Tab 2 — Komponen harmonik
        tab_comp = tk.Frame(nb, bg=BG)
        nb.add(tab_comp, text="  🎵  Komponen Harmonik  ")
        self._build_comp_tab(tab_comp)

        # Tab 3 — Tabel data
        tab_tbl = tk.Frame(nb, bg=BG)
        nb.add(tab_tbl, text="  📋  Tabel Data  ")
        self._build_table_tab(tab_tbl)

        # Tab Validasi
        tab_val = tk.Frame(nb, bg=BG)
        nb.add(tab_val, text="  ✅  Validasi  ")
        self.fig_val = Figure(figsize=(8, 5), facecolor=BG)
        self.canvas_val = FigureCanvasTkAgg(self.fig_val, tab_val)
        self.canvas_val.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.nb = nb
        self.tab_val = tab_val

        # Tab 4 — Cara Penggunaan
        tab_met = tk.Frame(nb, bg=BG)
        nb.add(tab_met, text="  📋  Cara Penggunaan  ")
        self._build_method_tab(tab_met)

        return frame

    def _build_plot_tab(self, parent):
        self.fig = Figure(figsize=(8, 4.5), facecolor=BG)
        self.ax  = self.fig.add_subplot(111, facecolor=BG2)
        self._style_ax(self.ax)
        self.ax.set_title("Muat file dan jalankan analisis untuk menampilkan grafik",
                          color=TXT2, fontsize=9)
        self.canvas_plot = FigureCanvasTkAgg(self.fig, parent)
        self.canvas_plot.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        tb = NavigationToolbar2Tk(self.canvas_plot, parent)
        tb.config(background=BG3)
        for w in tb.winfo_children():
            try: w.config(background=BG3, foreground=TXT)
            except: pass
        tb.update()

    def _build_comp_tab(self, parent):
        self.fig_comp = Figure(figsize=(8,4.5), facecolor=BG)
        self.canvas_comp = FigureCanvasTkAgg(self.fig_comp, parent)
        self.canvas_comp.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def _build_table_tab(self, parent):
        cols = ("No","Datetime","Water Level (m)","Keterangan")
        self.tree = ttk.Treeview(parent, columns=cols, show='headings',
                                 style='Tidal.Treeview')
        style = ttk.Style()
        style.configure('Tidal.Treeview', background=BG2, foreground=TXT,
                        fieldbackground=BG2, rowheight=22, font=("Consolas",9))
        style.configure('Tidal.Treeview.Heading', background=BG3,
                        foreground=ACCENT, font=("Segoe UI",9,"bold"))
        style.map('Tidal.Treeview', background=[('selected', ACCENT)],
                  foreground=[('selected', 'white')])
        for col in cols:
            w = 80 if col=="No" else (180 if col=="Datetime" else 140)
            self.tree.heading(col, text=col)
            self.tree.column(col, width=w, anchor=tk.CENTER)
        sb = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        self.tree.pack(fill=tk.BOTH, expand=True)

    def _build_method_tab(self, parent):
        txt = tk.Text(parent, bg=BG2, fg=TXT, font=("Segoe UI",9),
                      relief=tk.FLAT, wrap=tk.WORD, padx=16, pady=12,
                      insertbackground=TXT, state=tk.DISABLED)
        sb = ttk.Scrollbar(parent, orient=tk.VERTICAL, command=txt.yview)
        txt.configure(yscrollcommand=sb.set)
        sb.pack(side=tk.RIGHT, fill=tk.Y)
        txt.pack(fill=tk.BOTH, expand=True)
        txt.configure(state=tk.NORMAL)
        txt.insert(tk.END, METODOLOGI_TEXT)
        txt.configure(state=tk.DISABLED)
        self.method_text = txt

    # ══════════════════════════════════════════
    # LOGIC
    # ══════════════════════════════════════════
    # ── Koreksi datum DEMNAS ──
    def _calc_dem(self):
        try:
            msl = float(self.msl_portal.get())
            n_ina = float(self.n_inageoid.get())
            lat, lon = float(self.lat_site.get()), float(self.lon_site.get())
        except Exception:
            messagebox.showwarning("Peringatan", "Isi MSL portal, N InaGeoid2020, lat dan lon dengan angka.")
            return
        manual = self.n_egm2008.get().strip()
        if manual:
            try:
                self._finish_dem(msl, n_ina, float(manual), "isian manual")
            except ValueError:
                messagebox.showwarning("Peringatan", "N EGM2008 harus angka (atau kosongkan agar dicari otomatis).")
            return
        self.dem_info.config(text="Mencari N EGM2008...", fg=WARN)
        self._set_status("Mengambil undulasi EGM2008...", None)

        def work():
            try:
                n_egm, src = egm2008_undulation(lat, lon)
                self.after(0, lambda: (self.n_egm2008.set(f"{n_egm:.4f}"),
                                       self._finish_dem(msl, n_ina, n_egm, src)))
            except Exception as e:
                msg = str(e)
                self.after(0, lambda: self.dem_info.config(
                    text="Gagal mengambil N EGM2008 otomatis (cek internet).\n"
                         "Isi manual dari kalkulator geoid, lalu klik lagi.\n" + msg[:120],
                    fg=ERR))
        threading.Thread(target=work, daemon=True).start()

    def _finish_dem(self, msl, n_ina, n_egm, src):
        off = msl + (n_ina - n_egm)
        self.dem_offset.set(round(off, 4))
        self.dem_source = src
        self.dem_info.config(
            text=f"Koreksi = {msl:.4f} + ({n_ina:.3f} − {n_egm:.3f}) = {off:.3f} m",
            fg=ACCENT2)
        self._set_status(f"Koreksi ke DEMNAS = {off:.3f} m", None)

    # ── Waktu & puncak ──
    def _tz(self):
        """Selisih waktu lokal terhadap UTC."""
        try:
            return timedelta(hours=float(self.utc_offset.get()))
        except Exception:
            return timedelta(0)

    def _tz_label(self):
        try:
            h = float(self.utc_offset.get())
        except Exception:
            h = 0.0
        return f"UTC{'+' if h >= 0 else '-'}{abs(h):g}"

    def _to_local(self, t_utc):
        return t_utc + self._tz()

    def _find_peaks(self):
        """Puncak lokal (pasang tinggi) dari data BIG, urut z menurun.
        Hasil: [(z, t_utc), ...]. Dataran puncak (nilai sama) dihitung sekali."""
        zvals = [r[1] for r in self.raw_data]
        times = [r[0] for r in self.raw_data]
        peaks = []
        for i in range(1, len(zvals) - 1):
            if zvals[i] > zvals[i-1] and zvals[i] >= zvals[i+1]:
                peaks.append((zvals[i], times[i]))
        peaks.sort(key=lambda x: x[0], reverse=True)
        return peaks

    def _load_file(self):
        path = filedialog.askopenfilename(
            title="Pilih file prediksi pasut BIG SRGI",
            filetypes=[("Text/CSV files","*.txt *.csv"), ("All files","*.*")])
        if not path:
            return
        self._set_status("Membaca file...", 10)
        try:
            rows = []
            with open(path, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    m = re.match(
                        r'^\s*([\d.]+)\s+([\d.]+)\s+'
                        r'(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})\s+([-\d.]+)',
                        line)
                    if m:
                        dt = datetime.strptime(m.group(3), '%Y-%m-%d %H:%M:%S')
                        rows.append((dt, float(m.group(4))))
                        if len(rows) == 1:
                            first_latlon = (float(m.group(1)), float(m.group(2)))
            if len(rows) < 24:
                messagebox.showerror("Error",
                    "Data terlalu sedikit atau format tidak dikenali.\n"
                    "Pastikan file dari menu 'Unduh Prediksi Pasut' BIG SRGI.")
                return
            self.raw_data = rows
            self.lat_site.set(round(first_latlon[0], 6))
            self.lon_site.set(round(first_latlon[1], 6))
            self.n_egm2008.set("")
            fname = os.path.basename(path)
            self.file_label.config(text=fname, fg=ACCENT2)
            t0, t1 = rows[0][0], rows[-1][0]
            zvals = [r[1] for r in rows]
            l0, l1 = self._to_local(t0), self._to_local(t1)
            info = (f"✓  {len(rows)} baris  |  "
                    f"{l0.strftime('%Y-%m-%d')} – {l1.strftime('%Y-%m-%d')} ({self._tz_label()})\n"
                    f"z: min={min(zvals):.3f} m, max={max(zvals):.3f} m")
            self.data_info.config(text=info, fg=ACCENT2)
            self.start_dt_str.set(l0.strftime('%Y-%m-%d %H:%M'))
            self._set_status(f"File dimuat: {len(rows)} data.", 30)
        except Exception as e:
            messagebox.showerror("Error membaca file", str(e))
            self._set_status("Gagal membaca file.", 0)

    def _run_analysis(self):
        if not self.raw_data:
            messagebox.showwarning("Peringatan", "Muat file BIG SRGI terlebih dahulu.")
            return
        if not UTIDE_OK:
            messagebox.showerror("Error",
                "Library 'utide' tidak ditemukan.\n"
                "Install dengan: pip install utide")
            return
        self._set_status("Menjalankan analisis harmonik...", 40)
        threading.Thread(target=self._do_analysis, daemon=True).start()

    def _do_analysis(self):
        try:
            times_dt = [r[0] for r in self.raw_data]
            water    = np.array([r[1] for r in self.raw_data])
            t_np     = np.array([np.datetime64(t) for t in times_dt])

            coef = utide.solve(
                t_np, water,
                lat=self.lat_site.get(),
                constit='auto',          # semua konstituen yang terpisah (kriteria Rayleigh)
                method='ols',
                conf_int='linear',
                verbose=False)
            self.coef = coef
            self.n_const = len(coef.name)
            self.z0 = float(np.asarray(coef['mean']).ravel()[0])

            names = list(coef.name)
            T_dict = {'M2':12.4206,'S2':12.0000,'K1':23.9345,'O1':25.8194,
                      'N2':12.6583,'K2':11.9672,'P1':24.0659,'Q1':26.8684}
            self.constituents = {}
            for name in ['M2','S2','K1','O1','N2','K2','P1','Q1']:
                if name in names:
                    idx = names.index(name)
                    self.constituents[name] = {
                        'A': float(coef.A[idx]),
                        'g': float(coef.g[idx]),
                        'T': T_dict.get(name, 0)
                    }

            self.after(0, self._update_const_display)
            self.after(0, self._plot_comp_chart)
            self.after(0, lambda: self._set_status(
                f"Analisis harmonik selesai: {self.n_const} konstituen (otomatis). "
                f"Siap generate time series.", 60))
        except Exception as e:
            self.after(0, lambda: messagebox.showerror("Error analisis", str(e)))
            self.after(0, lambda: self._set_status("Analisis gagal.", 0))

    def _update_const_display(self):
        for w in self.const_frame.winfo_children():
            w.destroy()
        if not self.constituents:
            return

        main4 = ['M2','S2','K1','O1']
        tk.Label(self.const_frame,
                 text=f"{getattr(self, 'n_const', 0)} konstituen dipakai (otomatis); "
                      f"ditampilkan 4 utama",
                 font=("Segoe UI",8,"italic"), bg=BG2, fg=ACCENT2).pack(anchor=tk.W)
        header = tk.Frame(self.const_frame, bg=BG3)
        header.pack(fill=tk.X, pady=(2,1))
        for h, w in [("Komp",4),("A (m)",8),("Fase (°)",10),("T (jam)",9)]:
            tk.Label(header, text=h, font=("Segoe UI",8,"bold"),
                     bg=BG3, fg=ACCENT, width=w, anchor=tk.CENTER).pack(side=tk.LEFT)

        for name in main4:
            if name not in self.constituents:
                continue
            c = self.constituents[name]
            row = tk.Frame(self.const_frame, bg=BG2)
            row.pack(fill=tk.X, pady=1)
            clr = ACCENT if name in ['M2','S2'] else ACCENT2
            tk.Label(row, text=name, font=("Consolas",9,"bold"),
                     bg=BG2, fg=clr, width=4, anchor=tk.CENTER).pack(side=tk.LEFT)
            tk.Label(row, text=f"{c['A']:.4f}", font=("Consolas",9),
                     bg=BG2, fg=TXT, width=8, anchor=tk.CENTER).pack(side=tk.LEFT)
            tk.Label(row, text=f"{c['g']:.2f}°", font=("Consolas",9),
                     bg=BG2, fg=TXT, width=10, anchor=tk.CENTER).pack(side=tk.LEFT)
            tk.Label(row, text=f"{c['T']:.4f}", font=("Consolas",9),
                     bg=BG2, fg=TXT2, width=9, anchor=tk.CENTER).pack(side=tk.LEFT)

        # Form number
        if all(k in self.constituents for k in ['M2','S2','K1','O1']):
            F = ((self.constituents['K1']['A'] + self.constituents['O1']['A']) /
                 (self.constituents['M2']['A'] + self.constituents['S2']['A']))
            tipe = ("Semidiurnal" if F < 0.25 else
                    "Campuran dominan semidiurnal" if F < 1.5 else
                    "Campuran dominan diurnal" if F < 3.0 else "Diurnal")
            sep = tk.Frame(self.const_frame, bg=BORDER, height=1)
            sep.pack(fill=tk.X, pady=4)
            tk.Label(self.const_frame,
                     text=f"Form Number F = {F:.4f}  →  {tipe}",
                     font=("Segoe UI",8,"bold"), bg=BG2, fg=WARN,
                     wraplength=300, justify=tk.LEFT).pack(anchor=tk.W)

    def _generate(self):
        if self.coef is None:
            messagebox.showwarning("Peringatan",
                "Jalankan Analisis Harmonik terlebih dahulu.")
            return
        self._set_status("Generating time series...", 70)
        threading.Thread(target=self._do_generate, daemon=True).start()

    def _do_generate(self):
        try:
            hat = self.hat_val.get()
            lat = self.lat_val.get()
            dur = self.duration_h.get()
            ivl = self.interval_min.get()
            # Waktu mulai diisi dalam waktu lokal; utide bekerja dalam UTC
            start = datetime.strptime(self.start_dt_str.get().strip(),
                                      '%Y-%m-%d %H:%M')
            start_utc = start - self._tz()

            n_steps = int(dur * 60 / ivl) + 1
            t_out = np.array([
                np.datetime64(start_utc + timedelta(minutes=i*ivl))
                for i in range(n_steps)])

            tide_out = utide.reconstruct(t_out, self.coef, verbose=False)
            raw = np.asarray(tide_out.h, dtype=float)

            # Koreksi datum: eta_MSL = h - Z0 (Z0 = muka air rata-rata data BIG)
            scaled = [round(float(v - self.z0), 4) for v in raw]
            metode = f"koreksi datum (Z0 = {self.z0:.3f} m)"
            self.method_used = metode
            self.start_utc = start_utc
            self.end_utc   = start_utc + timedelta(minutes=(n_steps - 1) * ivl)
            self.val = None

            self.scaled_values = scaled
            self.time_labels = [
                (start + timedelta(minutes=i*ivl)).strftime('%Y-%m-%d %H:%M')
                for i in range(n_steps)]

            self.after(0, self._update_plot)
            self.after(0, self._update_table)
            self.after(0, lambda: self._set_status(
                f"Selesai: {n_steps} titik, interval {ivl} menit, durasi {dur} jam | "
                f"{metode} | maks={max(scaled):.3f} m, min={min(scaled):.3f} m",
                100))
        except Exception as e:
            self.after(0, lambda: messagebox.showerror("Error generate", str(e)))
            self.after(0, lambda: self._set_status("Generate gagal.", 0))

    # ══════════════════════════════════════════
    # VALIDASI
    # ══════════════════════════════════════════
    def _run_validation(self):
        if self.coef is None or not self.scaled_values:
            messagebox.showwarning("Peringatan",
                "Jalankan analisis harmonik dan generate time series terlebih dahulu.")
            return
        self._set_status("Menjalankan validasi...", 85)
        threading.Thread(target=self._do_validation, daemon=True).start()

    def _do_validation(self):
        try:
            times = [r[0] for r in self.raw_data]
            z     = [r[1] for r in self.raw_data]
            self.val = compute_validation(times, z, self.coef, self.z0,
                                          self.start_utc, self.end_utc)
            self.after(0, self._show_validation)
        except Exception as e:
            self.after(0, lambda: messagebox.showerror("Error validasi", str(e)))
            self.after(0, lambda: self._set_status("Validasi gagal.", 0))

    def _show_validation(self):
        f, w = self.val['full'], self.val['win']
        lines = ["Setahun penuh (n={}):".format(f['n']),
                 "  RMSE={:.3f} m  R={:.4f}".format(f['RMSE'], f['R']),
                 "  NSE={:.4f}  Bias={:+.3f} m".format(f['NSE'], f['Bias'])]
        if w:
            lines += ["Jendela simulasi (n={}):".format(w['n']),
                      "  RMSE={:.3f} m  R={:.4f}".format(w['RMSE'], w['R']),
                      "  ΔPasang={:+.3f} m  ΔSurut={:+.3f} m".format(w['dHigh'], w['dLow'])]
        else:
            lines.append("Jendela simulasi di luar rentang data BIG.")
        self.val_info.config(text="\n".join(lines))
        self.fig_val.clear()
        self._draw_validation(self.fig_val, light=False)
        self.canvas_val.draw()
        self.nb.select(self.tab_val)
        self._set_status(f"Validasi selesai: RMSE setahun {f['RMSE']:.3f} m, R {f['R']:.4f}", 100)

    def _draw_validation(self, fig, light=False):
        if light:
            bg, txt, txt2, grid, c_mod, c_obs = "#FFFFFF", "#333333", "#666666", "#E0E4EC", "#2A7FD4", "#D94040"
        else:
            bg, txt, txt2, grid, c_mod, c_obs = BG2, TXT, TXT2, BORDER, TIDE_CLR, WARN
        gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.7, 1],
                               wspace=0.28, left=0.08, right=0.97, top=0.88, bottom=0.16)
        tz = self._tz()

        # (a) jendela simulasi: kurva rekonstruksi + titik BIG
        ax1 = fig.add_subplot(gs[0]); ax1.set_facecolor(bg)
        t_mod = [datetime.strptime(s, '%Y-%m-%d %H:%M') for s in self.time_labels]
        ax1.plot(t_mod, self.scaled_values, color=c_mod, lw=1.8, label='Rekonstruksi (15 menit)')
        times = [r[0] for r in self.raw_data]
        m = self.val['mask']
        if m.any():
            t_obs = [t + tz for t, k in zip(times, m) if k]
            ax1.plot(t_obs, self.val['obs'][m], 'o', ms=4, color=c_obs,
                     label='Data BIG (per jam)', zorder=3)
        w = self.val['win']
        ttl = "(a) Jendela simulasi"
        if w:
            ttl += f"  RMSE = {w['RMSE']:.3f} m, R = {w['R']:.3f}"
        ax1.set_title(ttl, color=txt, fontsize=9)
        ax1.set_xlabel(f"Waktu ({self._tz_label()})", color=txt2, fontsize=9)
        ax1.set_ylabel("Elevasi muka air (m, MSL)", color=txt2, fontsize=9)
        import matplotlib.dates as mdates
        ax1.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m\n%H:%M'))
        ax1.legend(fontsize=8, loc='best', facecolor=bg, edgecolor=grid, labelcolor=txt)

        # (b) setahun: sebaran BIG vs rekonstruksi
        ax2 = fig.add_subplot(gs[1]); ax2.set_facecolor(bg)
        obs, mod = self.val['obs'], self.val['mod']
        ax2.plot(obs, mod, '.', ms=2, alpha=0.35, color=c_mod)
        lo, hi = min(obs.min(), mod.min()) - 0.1, max(obs.max(), mod.max()) + 0.1
        ax2.plot([lo, hi], [lo, hi], '--', color=c_obs, lw=1.2, label='Garis 1:1')
        ax2.set_xlim(lo, hi); ax2.set_ylim(lo, hi); ax2.set_aspect('equal')
        f = self.val['full']
        ax2.set_title(f"(b) Setahun  RMSE = {f['RMSE']:.3f} m, R = {f['R']:.3f}",
                      color=txt, fontsize=9)
        ax2.set_xlabel("Data BIG (m, MSL)", color=txt2, fontsize=9)
        ax2.set_ylabel("Rekonstruksi (m, MSL)", color=txt2, fontsize=9)
        ax2.legend(fontsize=8, loc='upper left', facecolor=bg, edgecolor=grid, labelcolor=txt)

        for ax in (ax1, ax2):
            ax.tick_params(colors=txt2, labelsize=8)
            for sp in ('top', 'right'):
                ax.spines[sp].set_visible(False)
            for sp in ('bottom', 'left'):
                ax.spines[sp].set_color(grid)
            ax.grid(True, color=grid, lw=0.5, alpha=0.7, linestyle='--')

    def _export_validation(self):
        if self.val is None:
            messagebox.showwarning("Peringatan", "Jalankan validasi terlebih dahulu.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv", filetypes=[("CSV files", "*.csv")],
            initialfile="validasi_pasut_BIG_SRGI.csv")
        if not path:
            return
        try:
            f, w = self.val['full'], self.val['win']
            tz = self._tz()
            times = [r[0] for r in self.raw_data]
            with open(path, 'w', newline='', encoding='utf-8-sig') as fh:
                wr = csv.writer(fh)
                wr.writerow(["# Validasi rekonstruksi utide terhadap data prediksi BIG SRGI"])
                wr.writerow([f"# Datum MSL (Z0 = {self.z0:.4f} m); zona waktu {self._tz_label()}"])
                wr.writerow(["# statistik", "periode", "n", "RMSE_m", "MAE_m", "Bias_m", "R", "NSE", "galat_maks_m"])
                for name, d in (("setahun", f), ("jendela", w)):
                    if d:
                        wr.writerow(["#", name, d['n'], f"{d['RMSE']:.4f}", f"{d['MAE']:.4f}",
                                     f"{d['Bias']:.4f}", f"{d['R']:.5f}", f"{d['NSE']:.5f}",
                                     f"{d['max_err']:.4f}"])
                if w:
                    wr.writerow([f"# Selisih pasang tertinggi = {w['dHigh']:+.4f} m; "
                                 f"surut terendah = {w['dLow']:+.4f} m"])
                wr.writerow(["datetime", "BIG_m_MSL", "rekonstruksi_m_MSL", "selisih_m"])
                for t, k, o, mo in zip(times, self.val['mask'], self.val['obs'], self.val['mod']):
                    if k:
                        wr.writerow([(t + tz).strftime('%Y-%m-%d %H:%M'),
                                     f"{o:.4f}", f"{mo:.4f}", f"{mo - o:+.4f}"])
            png = os.path.splitext(path)[0] + ".png"
            fig_exp = Figure(figsize=(12, 5), facecolor="#FFFFFF", dpi=300)
            self._draw_validation(fig_exp, light=True)
            fig_exp.savefig(png, dpi=300, bbox_inches='tight', facecolor="#FFFFFF")
            self.export_info.config(
                text=f"✓  Disimpan: {os.path.basename(path)} + {os.path.basename(png)}", fg=ACCENT2)
            self._set_status(f"Validasi diekspor: {path}", 100)
        except Exception as e:
            messagebox.showerror("Error export validasi", str(e))

    def _make_time_axis_labels(self):
        """Buat label x-axis: jam:menit saja, dengan penanda hari jika >1 hari."""
        labels = []
        start_dt = datetime.strptime(self.time_labels[0], '%Y-%m-%d %H:%M')
        for lbl in self.time_labels:
            dt = datetime.strptime(lbl, '%Y-%m-%d %H:%M')
            delta_days = (dt.date() - start_dt.date()).days
            time_str = dt.strftime('%H:%M')
            if delta_days > 0:
                labels.append(f"D+{delta_days}\n{time_str}")
            else:
                labels.append(time_str)
        return labels

    def _draw_tidal_chart(self, ax, light=False):
        """Gambar grafik pasut pada ax yang diberikan. light=True untuk background cerah."""
        n   = len(self.scaled_values)
        hat = self.hat_val.get()
        lat = self.lat_val.get()
        x   = list(range(n))

        if light:
            bg_plot  = "#FFFFFF"
            bg_fig   = "#F7F9FC"
            clr_tide = "#2A7FD4"
            clr_hat  = "#D94040"
            clr_lat  = "#D94040"
            clr_msl  = "#999999"
            clr_grid = "#E0E4EC"
            clr_txt  = "#333333"
            clr_txt2 = "#666666"
            clr_spn  = "#CCCCCC"
            alpha_fill = 0.18
        else:
            bg_plot  = BG2
            clr_tide = TIDE_CLR
            clr_hat  = HAT_CLR
            clr_lat  = LAT_CLR
            clr_msl  = MSL_CLR
            clr_grid = BORDER
            clr_txt  = TXT
            clr_txt2 = TXT2
            clr_spn  = BORDER
            alpha_fill = 0.15

        ax.set_facecolor(bg_plot)
        ax.fill_between(x, self.scaled_values, lat,
                        where=[v > lat for v in self.scaled_values],
                        alpha=alpha_fill, color=clr_tide)
        ax.plot(x, self.scaled_values, color=clr_tide,
                lw=2.0, label='Muka Air (m)', zorder=3)
        ax.axhline(hat, color=clr_hat, lw=1.3, ls='--', alpha=0.85,
                   label=f'HAT = {hat:.3f} m')
        ax.axhline(lat, color=clr_lat, lw=1.3, ls='--', alpha=0.85,
                   label=f'LAT = {lat:.3f} m')
        ax.axhline(0, color=clr_msl, lw=1.0, ls=':', alpha=0.7,
                   label='MSL = 0.000 m')

        # X-axis: jam:menit saja
        time_only = self._make_time_axis_labels()
        step = max(1, n // 14)
        ticks = list(range(0, n, step))
        ax.set_xticks(ticks)
        ax.set_xticklabels([time_only[i] for i in ticks],
                           fontsize=8, color=clr_txt2)

        ax.set_ylabel("Elevasi muka air (m, MSL)", color=clr_txt2, fontsize=9)
        ax.set_xlabel(f"Waktu ({self._tz_label()})", color=clr_txt2, fontsize=9)

        # Judul dengan tanggal di subtitle
        start_lbl = self.time_labels[0]
        end_lbl   = self.time_labels[-1]
        ivl = self.interval_min.get()
        dur = self.duration_h.get()
        ax.set_title(
            f"Time Series Pasut Sintetis — BIG SRGI\n"
            f"{start_lbl}  s/d  {end_lbl}   |   Durasi {dur} jam, Interval {ivl} menit",
            color=clr_txt, fontsize=9, pad=8, linespacing=1.5)

        ax.legend(fontsize=8,
                  facecolor="#F0F4FA" if light else BG3,
                  edgecolor=clr_spn,
                  labelcolor=clr_txt,
                  loc='upper right',
                  framealpha=0.9)
        ax.set_ylim(min(lat, min(self.scaled_values)) - 0.3,
                    max(hat, max(self.scaled_values)) + 0.3)
        ax.tick_params(colors=clr_txt2, labelsize=8)
        ax.spines['bottom'].set_color(clr_spn)
        ax.spines['left'].set_color(clr_spn)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, color=clr_grid, lw=0.6, alpha=0.7, linestyle='--')

    def _update_plot(self):
        self.ax.clear()
        self._style_ax(self.ax)
        if not self.scaled_values:
            return
        self._draw_tidal_chart(self.ax, light=False)
        self.canvas_plot.draw()

    def _plot_comp_chart(self):
        self.fig_comp.clear()
        gs = gridspec.GridSpec(2, 1, figure=self.fig_comp,
                               hspace=0.45, top=0.92, bottom=0.12)

        # Bar chart amplitudo
        ax1 = self.fig_comp.add_subplot(gs[0])
        ax1.set_facecolor(BG2)
        names  = list(self.constituents.keys())
        amps   = [self.constituents[n]['A'] for n in names]
        colors = [ACCENT if n in ['M2','S2'] else ACCENT2 for n in names]
        bars = ax1.bar(names, amps, color=colors, width=0.5, zorder=3)
        ax1.bar_label(bars, fmt='%.3f', fontsize=8, color=TXT, padding=2)
        self._style_ax(ax1)
        ax1.set_title(f"Amplitudo 8 Komponen Utama (m) — dari {getattr(self, 'n_const', 0)} konstituen",
                      color=TXT, fontsize=9)
        ax1.set_ylabel("Amplitudo (m)", color=TXT2, fontsize=8)

        # Fase
        ax2 = self.fig_comp.add_subplot(gs[1])
        ax2.set_facecolor(BG2)
        phases = [self.constituents[n]['g'] for n in names]
        bars2  = ax2.bar(names, phases, color=WARN, width=0.5, zorder=3)
        ax2.bar_label(bars2, fmt='%.1f°', fontsize=8, color=TXT, padding=2)
        self._style_ax(ax2)
        ax2.set_title("Fase Komponen Harmonik (°)", color=TXT, fontsize=9)
        ax2.set_ylabel("Fase (°)", color=TXT2, fontsize=8)
        ax2.set_ylim(0, 400)

        self.canvas_comp.draw()

    def _update_table(self):
        for row in self.tree.get_children():
            self.tree.delete(row)
        hat = self.hat_val.get()
        lat = self.lat_val.get()
        for i, (dt_str, val) in enumerate(zip(self.time_labels, self.scaled_values), 1):
            ket = ""
            if abs(val - hat) < 0.01:   ket = "≈ HAT"
            elif abs(val - lat) < 0.01: ket = "≈ LAT"
            elif abs(val) < 0.05:       ket = "≈ MSL"
            elif val > hat * 0.9:       ket = "Pasang tinggi"
            elif val < lat * 0.9:       ket = "Surut rendah"
            tag = 'high' if val > 0.5 else ('low' if val < -0.5 else 'mid')
            self.tree.insert('', tk.END,
                             values=(i, dt_str, f"{val:.4f}", ket),
                             tags=(tag,))
        self.tree.tag_configure('high', foreground=ACCENT)
        self.tree.tag_configure('low',  foreground=WARN)
        self.tree.tag_configure('mid',  foreground=TXT)


    def _find_spring_tide(self):
        if not self.raw_data:
            messagebox.showwarning("Peringatan", "Muat file BIG SRGI terlebih dahulu.")
            return
        try:
            rl_hours = float(self.rising_limb_h.get())
        except:
            messagebox.showwarning("Peringatan", "Masukkan nilai rising limb BC hulu yang valid (jam).")
            return

        peaks = self._find_peaks()
        if not peaks:
            self.spring_result.config(text="Tidak ditemukan puncak pasut.")
            return

        top5 = [(z, self._to_local(t)) for z, t in peaks[:5]]
        best_time = top5[0][1]
        sim_start = best_time - timedelta(hours=rl_hours)
        self.start_dt_str.set(sim_start.strftime('%Y-%m-%d %H:%M'))

        lines = [f"Rising limb BC hulu: {rl_hours:.1f} jam",
                 f"Top 5 Spring Tide (z tertinggi, {self._tz_label()}):"]
        for i, (z, t) in enumerate(top5, 1):
            lines.append(f"{i}. {t.strftime('%Y-%m-%d %H:%M')}  z={z:.3f} m")
        lines.append(f"\n✓ Waktu mulai simulasi:")
        lines.append(f"  {sim_start.strftime('%Y-%m-%d %H:%M')}")
        lines.append(f"  (= puncak pasut − {rl_hours:.1f} jam rising limb)")
        self.spring_result.config(text="\n".join(lines), fg=WARN)
        self._set_status(
            f"Spring tide: {best_time.strftime('%Y-%m-%d %H:%M')} | "
            f"Mulai simulasi: {sim_start.strftime('%Y-%m-%d %H:%M')}", 80)

    def _show_top_z(self):
        if not self.raw_data:
            messagebox.showwarning("Peringatan", "Muat file BIG SRGI terlebih dahulu.")
            return
        # Hanya puncak lokal, agar satu pasang tinggi tidak muncul berulang
        top10 = [(z, self._to_local(t)) for z, t in self._find_peaks()[:10]]

        # Buat popup window
        win = tk.Toplevel(self)
        win.title("Top 10 Muka Air Tertinggi — BIG SRGI")
        win.configure(bg=BG)
        win.geometry("520x420")
        win.resizable(False, False)

        tk.Label(win, text=f"Top 10 Puncak Pasang Tertinggi ({self._tz_label()})",
                 font=("Segoe UI", 11, "bold"), bg=BG, fg=TXT).pack(pady=(16,4))
        rl_hours = 6.0
        try:
            rl_hours = float(self.rising_limb_h.get())
        except:
            pass
        tk.Label(win,
                 text=f"Klik baris untuk set waktu mulai simulasi  (puncak − {rl_hours:.1f} jam rising limb BC hulu)",
                 font=("Segoe UI", 8), bg=BG, fg=TXT2).pack(pady=(0,10))

        frame_tbl = tk.Frame(win, bg=BG)
        frame_tbl.pack(fill=tk.BOTH, expand=True, padx=16)

        cols = ("No", "Tanggal & Waktu", "z (m)", "Gunakan sebagai")
        tree = ttk.Treeview(frame_tbl, columns=cols, show='headings', height=10)
        style = ttk.Style()
        style.configure('Top.Treeview', background=BG2, foreground=TXT,
                        fieldbackground=BG2, rowheight=28, font=("Consolas", 9))
        style.configure('Top.Treeview.Heading', background=BG3,
                        foreground=ACCENT, font=("Segoe UI", 9, "bold"))
        style.map('Top.Treeview', background=[('selected', ACCENT)],
                  foreground=[('selected', 'white')])
        tree.configure(style='Top.Treeview')

        widths = [40, 180, 90, 160]
        for col, w in zip(cols, widths):
            tree.heading(col, text=col)
            tree.column(col, width=w, anchor=tk.CENTER)

        for i, (z, t) in enumerate(top10, 1):
            sim_start = t - timedelta(hours=rl_hours)
            tree.insert('', tk.END,
                        values=(i, t.strftime('%Y-%m-%d %H:%M'),
                                f"{z:.4f}", sim_start.strftime('%Y-%m-%d %H:%M')),
                        tags=(str(i),))

        tree.tag_configure('1', foreground="#B8860B")
        tree.tag_configure('2', foreground="#6B7280")
        tree.tag_configure('3', foreground="#A0522D")

        def on_select(event):
            sel = tree.selection()
            if not sel:
                return
            vals = tree.item(sel[0], 'values')
            sim_start_str = vals[3]
            self.start_dt_str.set(sim_start_str)
            info_lbl.config(
                text=f"✓ Waktu mulai diset: {sim_start_str}",
                fg=ACCENT2)
            self.spring_result.config(
                text=f"Waktu mulai: {sim_start_str}\n"
                     f"(puncak z={vals[2]} m pada {vals[1]} − {rl_hours:.1f} jam rising limb)",
                fg=WARN)

        tree.bind('<<TreeviewSelect>>', on_select)
        tree.pack(fill=tk.BOTH, expand=True)

        info_lbl = tk.Label(win, text="", font=("Segoe UI", 9, "bold"),
                            bg=BG, fg=ACCENT2)
        info_lbl.pack(pady=8)

        self._btn(win, "Tutup", win.destroy, color=BTN_NEUTRAL).pack(pady=(0,12))

    def _export_csv(self):
        if not self.scaled_values:
            messagebox.showwarning("Peringatan", "Generate time series terlebih dahulu.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files","*.csv")],
            initialfile="tidal_timeseries_BIG_SRGI.csv")
        if not path:
            return
        try:
            with open(path, 'w', newline='', encoding='utf-8-sig') as f:
                w = csv.writer(f)
                w.writerow(["# Synthetic Tidal Time Series — BIG SRGI"])
                w.writerow([f"# Lokasi: Lat {self.lat_site.get():.4f}°"])
                w.writerow([f"# HAT={self.hat_val.get()}, LAT={self.lat_val.get()}, MSL=0.000"])
                w.writerow([f"# Datum: MSL; metode: {getattr(self, 'method_used', '-')}"])
                w.writerow([f"# Zona waktu kolom datetime: {self._tz_label()} (data BIG asli dalam UTC)"])
                w.writerow([f"# Koreksi datum: MSL portal={self.msl_portal.get()}, N InaGeoid2020={self.n_inageoid.get()}, "
                            f"N EGM2008={self.n_egm2008.get() or '-'} ({getattr(self, 'dem_source', '-')})"])
                w.writerow([f"# Tinggi MSL terhadap datum DEM = {self.dem_offset.get():.4f} m "
                            f"(kolom water_level_m_DEM = water_level_m_MSL + nilai ini; pakai kolom ini di HEC-RAS)"])
                w.writerow([f"# Analisis harmonik utide (Codiga 2011): {getattr(self, 'n_const', 0)} konstituen (otomatis)"])
                w.writerow([f"# Komponen utama:"])
                for name, c in self.constituents.items():
                    w.writerow([f"#   {name}: A={c['A']:.4f} m, fase={c['g']:.2f}°"])
                off = float(self.dem_offset.get())
                w.writerow(["datetime", "water_level_m_MSL", "water_level_m_DEM"])
                for dt_str, val in zip(self.time_labels, self.scaled_values):
                    w.writerow([dt_str, val, round(val + off, 4)])
            self.export_info.config(
                text=f"✓  Disimpan: {os.path.basename(path)}", fg=ACCENT2)
            self._set_status(f"CSV disimpan: {path}", 100)
        except Exception as e:
            messagebox.showerror("Error export", str(e))

    def _export_png(self):
        if not self.scaled_values:
            messagebox.showwarning("Peringatan", "Generate time series terlebih dahulu.")
            return
        path = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG image","*.png")],
            initialfile="tidal_timeseries_BIG_SRGI.png")
        if not path:
            return
        try:
            # Buat figure baru dengan background cerah untuk export
            fig_exp = Figure(figsize=(12, 5), facecolor="#F7F9FC", dpi=300)
            ax_exp  = fig_exp.add_subplot(111)
            self._draw_tidal_chart(ax_exp, light=True)
            fig_exp.tight_layout(pad=1.5)
            fig_exp.savefig(path, dpi=300, bbox_inches='tight',
                            facecolor="#F7F9FC")
            plt.close('all')
            self.export_info.config(
                text=f"✓  Grafik disimpan: {os.path.basename(path)}", fg=ACCENT2)
            self._set_status(f"Grafik (background cerah) disimpan: {path}", 100)
        except Exception as e:
            messagebox.showerror("Error export PNG", str(e))

    # ── Helpers ──
    def _style_ax(self, ax):
        ax.set_facecolor(BG2)
        ax.tick_params(colors=TXT2, labelsize=8)
        ax.spines['bottom'].set_color(BORDER)
        ax.spines['left'].set_color(BORDER)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True, color=BORDER, lw=0.5, alpha=0.5, linestyle='--')
        ax.yaxis.label.set_color(TXT2)
        ax.xaxis.label.set_color(TXT2)

    def _set_status(self, msg, pct=None):
        self.status_var.set(msg)
        if pct is not None:
            self.progress_var.set(pct)
        self.update_idletasks()


# ─────────────────────────────────────────────
# TEKS METODOLOGI
# ─────────────────────────────────────────────
METODOLOGI_TEXT = """
CARA PENGGUNAAN — TIDAL TIME SERIES GENERATOR (BIG SRGI)
=========================================================

━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LANGKAH A — UNDUH DATA PREDIKSI PASUT DARI BIG SRGI
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

1. Buka browser → kunjungi:
   https://srgi.big.go.id/tides/bbng

2. Pada peta interaktif, klik lokasi titik yang diinginkan
   (di laut, dekat pantai area studi Anda).
   → Pastikan koordinat jatuh di laut, bukan di daratan.
   → Jika muncul "Site is out of model grid OR land",
     geser koordinat sedikit ke arah laut (offshore).

3. Pada panel yang muncul, pilih menu:
   "Unduh Prediksi Pasut"

4. Isi formulir:
   - Tanggal mulai  : 1 Januari tahun yang diinginkan
   - Tanggal akhir  : 31 Desember tahun yang sama
     (minimal 1 tahun = ±8.760 data per jam)
   - Koordinat      : sesuai titik yang dipilih di peta

5. Klik tombol Download → simpan file .txt

   Format isi file hasil unduhan:
   Lat  Lon  yyyy-mm-dd hh:mm:ss (UTC)  z(m)
   (setiap baris = satu data elevasi muka air per jam)

CATATAN DATUM:
   Nilai z pada file prediksi BIG SRGI sudah mengacu ke MSL
   (muka air rata-rata Z0 hasil analisis ≈ 0). Sebagai
   pengaman, aplikasi tetap mengurangkan Z0:
       eta_MSL = h - Z0
   HAT/LAT hanya dipakai sebagai garis pembanding di grafik.

CATATAN ZONA WAKTU:
   Waktu pada file BIG SRGI adalah UTC. Aplikasi menampilkan
   dan menerima waktu dalam zona lokal sesuai isian
   "Zona waktu" di panel 2 (Ternate/WIT = UTC+9).


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LANGKAH B — MUAT FILE DAN JALANKAN ANALISIS HARMONIK
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Panel Kiri → 1. Load File
  • Klik "📂 Muat File BIG SRGI"
  • Pilih file .txt hasil unduhan dari SRGI BIG
  • Informasi data (jumlah baris, rentang tanggal,
    z min/maks) akan tampil di bawah tombol

Panel Kiri → 2. Parameter Lokasi & Datum
  • Latitude dan Longitude terisi otomatis dari file BIG
  • Isi HAT dan LAT dari portal BIG SRGI (Model Pasut,
    klik titik yang sama pada peta)
  • Zona waktu: WIT (Ternate) = 9, WITA = 8, WIB = 7
  • Lanjutkan dengan koreksi datum (Langkah B2)


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LANGKAH B2 — KOREKSI DATUM KE DEMNAS (PANEL 2)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

MENGAPA PERLU:
  Hasil aplikasi mengacu ke MSL (MSL = 0). DEMNAS mengacu ke
  geoid EGM2008. Portal pasut BIG memberi tinggi MSL terhadap
  geoid nasional InaGeoid2020. Tanpa koreksi, kondisi batas
  di HEC-RAS bisa terlalu rendah sekitar 2 m.

  Koreksi = MSL portal + (N InaGeoid2020 − N EGM2008)

LANGKAH:
  1. Buka portal SRGI (https://srgi.big.go.id), menu
     "Model Pasut" → klik titik yang SAMA dengan file BIG
     (atau Single Coordinate dengan lat/lon di panel 2).
     Catat nilai MSL.
       → Ketik di kolom "MSL portal SRGI (m)"

  2. Di portal yang sama, menu "Model Geoid" →
     Single Coordinate → isi lat/lon yang sama,
     Ellipsoid height = 0 → Calculate.
     Catat nilai Geoid Undulation (N).
       → Ketik di kolom "N InaGeoid2020 (m)"

  3. Biarkan kolom "N EGM2008 (m)" KOSONG, lalu klik
     "Hitung Koreksi ke DEMNAS".
       → Aplikasi mencari N EGM2008 sendiri (perlu internet)
       → Kolom "Koreksi ke DEM (m)" terisi otomatis
       → Rumus dan hasilnya tampil di bawah tombol

  4. Jika muncul pesan gagal (tidak ada internet):
     cari N EGM2008 di kalkulator geoid, misalnya
     https://geographiclib.sourceforge.io/cgi-bin/GeoidEval
     (masukkan lat lon), ketik di kolom "N EGM2008 (m)",
     lalu klik tombol lagi.

  CONTOH TERNATE (0.820547, 127.395267):
     MSL portal      = 1.3588 m
     N InaGeoid2020  = 72.397 m
     N EGM2008       = 71.784 m (otomatis)
     Koreksi         = 1.3588 + (72.397 − 71.784) = 1.972 m

  HASIL:
     CSV ekspor berisi dua kolom muka air:
       water_level_m_MSL  → untuk grafik/laporan
       water_level_m_DEM  → untuk HEC-RAS (= MSL + koreksi)


Panel Kiri → 3. Analisis Harmonik
  • Klik "▶ Jalankan Analisis Harmonik"
  • Proses memakan waktu beberapa detik
  • Konstituen dipilih otomatis oleh utide: semua yang dapat
    dipisahkan dari panjang data (kriteria Rayleigh). Data
    setahun per jam menghasilkan puluhan konstituen.
  • Yang ditampilkan: amplitudo & fase 8 komponen utama
    (M2, S2, K1, O1, N2, K2, P1, Q1) dan Form Number


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LANGKAH C — TENTUKAN WAKTU SIMULASI (SKENARIO)
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Panel Kiri → 4. Analisis Waktu Simulasi
  Digunakan untuk menentukan waktu mulai simulasi HEC-RAS
  berdasarkan kondisi pasut paling kritis.

  INPUT: Rising limb BC hulu (jam)
    → Durasi naik hidrograf dari HEC-HMS (jam 0 sampai
      puncak debit). Contoh: 11 jam untuk Q100.
    → Nilai ini diperoleh dari grafik hidrograf HEC-HMS.

  ┌─────────────────────────────────────────────────────┐
  │  PILIH SALAH SATU tombol sesuai skenario simulasi:  │
  └─────────────────────────────────────────────────────┘

  TOMBOL "🌕 Cari Spring Tide (z tertinggi)"
  → Gunakan untuk SKENARIO 2 atau SKENARIO 3
  → Mencari 1 waktu puncak pasut tertinggi (spring tide)
    dari seluruh data setahun secara otomatis.
  → Waktu mulai = t_puncak − rising limb → terisi di panel 5.

    Skenario 2 (Spring Tide):
      BC hilir = time series pasut mulai dari waktu ini.
      BC hulu  = hidrograf Q rencana, mulai waktu yang sama.
      Puncak pasut terjadi saat simulasi berjalan.

    Skenario 3 (Worst Case — paling kritis):
      Sama seperti Skenario 2, TETAPI di HEC-HMS, waktu
      mulai hidrograf diset agar puncak debit Q tiba
      di muara TEPAT SAAT pasut sedang di puncaknya.
      → Ini adalah kondisi banjir terbesar + pasut tertinggi
        terjadi bersamaan.

  ─────────────────────────────────────────────────────

  TOMBOL "📅 Tampilkan Top 10 z Tertinggi"
  → Gunakan untuk SKENARIO 1 (Validasi)
  → Menampilkan 10 puncak pasang tertinggi (satu baris
    per kejadian pasang, dalam waktu lokal).
  → Pilih tanggal yang mendekati tanggal kejadian banjir
    aktual di lapangan → klik baris → waktu mulai terisi.

  ─────────────────────────────────────────────────────

  ⚠ PERHATIAN:
    Kedua tombol hanya mengisi "Waktu mulai" di panel 5.
    Jika keduanya diklik, yang berlaku adalah klik TERAKHIR.
    Pilih SATU tombol saja sesuai skenario yang dijalankan,
    lalu langsung lanjut ke Generate Time Series.


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LANGKAH D — GENERATE TIME SERIES
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Panel Kiri → 5. Parameter Generate
  • Waktu mulai : terisi otomatis dari Langkah C,
                  atau diisi manual (format: YYYY-MM-DD HH:MM)
  • Durasi (jam): panjang simulasi, contoh: 30 jam
  • Interval    : resolusi waktu, contoh: 15 menit

  Klik "🌊 Generate Time Series"
  → Grafik time series pasut muncul di panel kanan
  → Tab "📈 Grafik Time Series" menampilkan kurva pasut
    dengan garis referensi HAT, MSL, LAT


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
LANGKAH E — EXPORT HASIL
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Panel Kiri → 6. Validasi
  • Klik "▶ Jalankan Validasi" setelah Generate
  • Rekonstruksi utide dibandingkan dengan data per jam BIG:
      - setahun penuh  → RMSE, MAE, bias, R, NSE
      - jendela simulasi → RMSE, R, selisih pasang tertinggi
        dan surut terendah
  • Tab "✅ Validasi" menampilkan (a) kurva + titik BIG pada
    jendela simulasi dan (b) diagram sebar setahun (garis 1:1)
  • Catatan: data pembanding adalah prediksi BIG sendiri,
    jadi validasi ini mengukur galat rekonstruksi harmonik
    dan perubahan interval, bukan galat terhadap observasi
    lapangan.

Panel Kiri → 7. Export
  • "Export CSV"          → file .csv siap diimpor ke HEC-RAS
                            sebagai Stage Hydrograph (BC hilir)
  • "Export Grafik (PNG)" → gambar grafik latar terang,
                            300 dpi, siap untuk laporan/jurnal
  • "Export Validasi"     → CSV (statistik + tabel per jam)
                            dan PNG 300 dpi dengan nama sama


━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
CATATAN PENGGUNAAN DI HEC-RAS 2D
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

  • Tipe BC hilir : Stage Hydrograph
  • Datum elevasi : gunakan kolom water_level_m_DEM
  • MSL TIDAK sama dengan datum DEMNAS. Lakukan koreksi
    datum di panel 2 (lihat Langkah B2) sebelum export.
  • BC hulu dan BC hilir harus memiliki rentang waktu
    yang sama persis (HEC-RAS tidak mengizinkan beda range)
  • Waktu mulai simulasi = waktu puncak pasut − rising limb
    BC hulu (hidrograf dari HEC-HMS)
"""


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────
if __name__ == "__main__":
    if not UTIDE_OK:
        root = tk.Tk()
        root.withdraw()
        messagebox.showerror(
            "Dependensi tidak ditemukan",
            "Library 'utide' tidak terinstall.\n\n"
            "Install dengan perintah:\n"
            "  pip install utide numpy pandas matplotlib\n\n"
            "Kemudian jalankan ulang aplikasi ini.")
        root.destroy()
    else:
        app = TidalApp()
        app.mainloop()
