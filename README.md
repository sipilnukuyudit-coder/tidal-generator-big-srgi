# Tidal Time Series Generator — BIG SRGI

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23296857.svg)](https://doi.org/10.5281/zenodo.23296857)

Aplikasi Python berantarmuka grafis untuk menyusun **kondisi batas hilir pasang surut** bagi pemodelan banjir pesisir HEC-RAS 2D, berbasis data prediksi pasang surut Badan Informasi Geospasial (portal SRGI).

*A Python GUI tool that converts BIG SRGI tide predictions into datum-consistent downstream boundary conditions for HEC-RAS 2D coastal flood models.*

## Fitur

- Membaca file **Unduh Prediksi Pasut** dari portal SRGI BIG (data per jam, UTC).
- Analisis harmonik ulang dengan **UTide** (pemilihan konstituen otomatis, koreksi nodal) dan klasifikasi tipe pasang surut (bilangan Formzahl).
- Penentuan **waktu mulai simulasi**: puncak pasang tertinggi dikurangi rising limb hidrograf hulu (HEC-HMS).
- Rekonstruksi deret waktu pada interval pilihan (bawaan 15 menit), ditampilkan dalam waktu lokal (WIT/WITA/WIB).
- **Koreksi datum ke DEMNAS (EGM2008)**: cukup masukkan MSL dan N InaGeoid2020 dari portal SRGI; N EGM2008 diambil otomatis.
- **Validasi** terhadap data BIG (RMSE, MAE, bias, R, NSE) untuk setahun penuh dan jendela simulasi.
- Ekspor CSV (kolom MSL dan DEMNAS) siap dipakai sebagai *Stage Hydrograph* HEC-RAS, serta grafik PNG 300 dpi.

## Instalasi

Membutuhkan Python 3.9 atau lebih baru.

```bash
pip install -r requirements.txt
python tidal_generator_BIG_SRGI.py
```

`pyproj` bersifat opsional. Tanpa `pyproj`, nilai N EGM2008 diambil dari kalkulator GeoidEval GeographicLib (perlu internet). Jika tidak ada internet, nilainya dapat diisi manual.

## Cara pakai singkat

1. Unduh prediksi pasang surut setahun dari https://srgi.big.go.id (menu Model Pasut → Unduh Prediksi Pasut).
2. Muat file tersebut di panel 1.
3. Di panel 2, isi HAT, LAT, MSL portal, dan N InaGeoid2020 (menu Model Geoid, *ellipsoid height* = 0), lalu klik **Hitung Koreksi ke DEMNAS**.
4. Jalankan analisis harmonik (panel 3), tentukan waktu simulasi (panel 4), lalu **Generate** (panel 5).
5. Jalankan validasi (panel 6) dan ekspor hasil (panel 7).

Panduan lengkap ada di tab **Cara Penggunaan** di dalam aplikasi.

## Data

Repositori ini tidak menyertakan data BIG. Silakan unduh langsung dari portal SRGI dan patuhi ketentuan penggunaan data BIG.

## Sitasi

Jika memakai aplikasi ini, mohon sitasi:

> Priambodo, Y. A. (2026). *Tidal Time Series Generator — BIG SRGI* (versi 1.0.0) [perangkat lunak]. https://doi.org/10.5281/zenodo.23296857

Serta makalah terkait:

> Priambodo, Y. A. (2026). Pengembangan aplikasi Python pembangkit deret waktu pasang surut dari data prediksi BIG SRGI untuk kondisi batas model banjir pesisir. [Nama jurnal].

## Lisensi

MIT License. Lihat file `LICENSE`.

## Kontak

Yudit Agus Priambodo — Program Studi Teknik Sipil, Universitas Muhammadiyah Maluku Utara (sipilummu.yudit@gmail.com)
