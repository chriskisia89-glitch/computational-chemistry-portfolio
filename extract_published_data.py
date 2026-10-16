"""Reproduce extraction of approximate image coordinates from the source SI.

Input is the separately downloaded PDF, not included in the package. Extraction
is locked to its checksum; the graphical data are not raw instrument records.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import sys
import numpy as np
from scipy.ndimage import distance_transform_edt, maximum_filter, label

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_MD5 = '448ef6d3fb603024d3abef4c65af75e4'


def extract(pdf_path):
    import fitz
    from PIL import Image
    import io
    blob = Path(pdf_path).read_bytes()
    if hashlib.md5(blob).hexdigest() != EXPECTED_MD5:
        raise ValueError('source PDF checksum differs from the audited version')
    doc = fitz.open(stream=blob, filetype='pdf')
    page = doc[4]  # S-5, Table S.2, C-PBDAZ-300.
    image = doc.extract_image(page.get_images()[0][0])['image']
    im = np.asarray(Image.open(io.BytesIO(image)).convert('RGB'))
    if im.shape != (1259, 1952, 3):
        raise ValueError('unexpected source raster size')
    black = im.max(axis=2) < 80

    def centers(values, offset):
        idx = np.flatnonzero(values)
        groups = np.split(idx, np.where(np.diff(idx)>1)[0]+1)
        return np.array([group.mean()+offset for group in groups])

    x_ticks = centers(black[582, 490:1180], 490)
    y_ticks = centers(black[25:580, 483], 25)
    if len(x_ticks) != 4 or len(y_ticks) != 7:
        raise ValueError('could not identify the expected axis ticks')
    x_affine = np.polyfit([0, 1, 2, 3], x_ticks, 1)
    y_affine = np.polyfit(np.arange(7, 3.5, -.5), y_ticks, 1)
    records = []
    for color, channel, temperature in [('blue', 2, 273), ('red', 0, 298)]:
        others = [i for i in range(3) if i != channel]
        mask = (im[:,:,channel]>160) & (im[:,:,others[0]]<100) & (im[:,:,others[1]]<100)
        mask[:30] = False; mask[580:] = False; mask[:,:480] = False; mask[:,1180:] = False
        dist = distance_transform_edt(mask)
        peaks = (dist == maximum_filter(dist, size=9)) & (dist > 5)
        labels, count = label(peaks)
        candidates=[]
        for k in range(1, count+1):
            yy, xx = np.where(labels == k)
            candidates.append((float(xx.mean()), float(yy.mean()), float(dist[yy,xx].max())))
        kept=[]
        for x, y, radius in sorted(candidates, key=lambda a: (-a[2], a[0], a[1])):
            # Dense red markers merge. Thin them to avoid treating every pixel
            # plateau as a separate observation. Preserve all blue candidates.
            if color == 'blue' or all(np.hypot(x-k[0], y-k[1]) >= 15 for k in kept):
                kept.append((x, y, radius))
        for x, y, radius in sorted(kept):
            records.append({'temperature_K': temperature,
                            'uptake_mmol_g': (x-x_affine[1])/x_affine[0],
                            'log_p_over_n_source': (y-y_affine[1])/y_affine[0],
                            'source_x_pixel': x, 'source_y_pixel': y,
                            'kind': 'digitized_plot_coordinate'})
    dest = ROOT/'data/published'; dest.mkdir(parents=True, exist_ok=True)
    with (dest/'c_pbdaz_300_digitized.csv').open('w', newline='', encoding='utf-8') as f:
        writer=csv.DictWriter(f, fieldnames=list(records[0])); writer.writeheader(); writer.writerows(records)
    provenance = {
      'paper_doi': '10.1021/acs.chemmater.3c01717',
      'supporting_information_doi': '10.1021/acs.chemmater.3c01717.s001',
      'figshare_article_id': 25019778, 'figshare_file_id': 44106648,
      'source_url': 'https://api.figshare.com/v2/articles/25019778',
      'download_url': 'https://ndownloader.figshare.com/files/44106648',
      'source_pdf_md5': EXPECTED_MD5, 'source_pdf_sha256': hashlib.sha256(blob).hexdigest(),
      'source_page': 'S-5', 'source_panel': 'Table S.2, C-PBDAZ-300, upper-left plot',
      'source_authors': ['Vaibhav A. Edlabadkar', 'Rushi U. Soni', 'A. B. M. Shaheen ud Doulah', 'Stephen Y. Owusu', 'Samuel Hackett', 'Joshua M. Bartels', 'Nicholas Leventis', 'Chariklia Sotiriou-Leventis'],
      'source_license': 'CC BY-NC 4.0', 'source_license_url': 'https://creativecommons.org/licenses/by-nc/4.0/',
      'raster_size_px': [1952,1259], 'x_tick_pixels': x_ticks.tolist(),
      'x_tick_values_mmol_g': [0,1,2,3], 'y_tick_pixels': y_ticks.tolist(),
      'y_tick_values': list(np.arange(7,3.5,-.5)),
      'x_pixel_equals_slope_times_loading_plus_intercept': x_affine.tolist(),
      'y_pixel_equals_slope_times_logordinate_plus_intercept': y_affine.tolist(),
      'point_count': len(records), 'counts_by_temperature': {str(t):sum(r['temperature_K']==t for r in records) for t in (273,298)},
      'data_type': 'Approximate coordinates extracted from a published raster plot, including markers and connected curves; not original instrument measurements.',
      'pressure_reference': 'The panel does not state its pressure reference unit. The source log ordinate is fitted directly; inferred absolute pressure is not labelled Torr or bar.',
      'digitization_limit': 'No experimental uncertainty is assigned. Finite pixel resolution, merged markers and plotted fit lines limit interpretation; no exact reproduction of the authors\' fit or reported errors is claimed.',
      'source_retrieved_utc': '2026-10-16',
      'algorithm': 'Threshold RGB channels; Euclidean distance transform; 9x9 local maxima above radius 5 px; centroid connected maxima; red candidates thinned to >=15 px separation.'}
    (dest/'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (dest/'published_coefficients.json').write_text(json.dumps({
      'sample':'C-PBDAZ-300','source':'Supporting Information S-5, Table S.2',
      'm':3,'n':1,'a':[-2466.79553,626.49851,-42.58868,2.78892],
      'b':[11.82664,-1.26531], 'reported_q0_kj_mol_rounded':21,
      'reported_r2':.99901,'reported_reduced_chi2':5.00934e-4},indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax=plt.subplots(figsize=(9,6))
    ax.imshow(im); ax.set_xlim(280,1210); ax.set_ylim(730,0)
    for temp, color in [(273,'cyan'),(298,'lime')]:
        rows=[r for r in records if r['temperature_K']==temp]
        ax.scatter([r['source_x_pixel'] for r in rows],[r['source_y_pixel'] for r in rows],facecolors='none',edgecolors=color,s=55,linewidths=1.0,label=f'{temp} K extracted coordinates')
    ax.set_xlabel('Source image x (pixels)');ax.set_ylabel('Source image y (pixels)')
    ax.set_title('Digitization audit: C-PBDAZ-300, Supporting Information S-5')
    ax.legend(loc='upper right',fontsize=8)
    fig.text(.01,.015,'Adapted from Edlabadkar et al., DOI 10.1021/acs.chemmater.3c01717.s001; source CC BY-NC 4.0.',fontsize=7)
    fig.tight_layout(rect=(0,.04,1,1))
    results=ROOT/'results';results.mkdir(exist_ok=True)
    fig.savefig(results/'digitization_audit.png',dpi=170);plt.close(fig)
    print(json.dumps({'rows':len(records),'counts':provenance['counts_by_temperature']},indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('source_pdf')
    extract(parser.parse_args().source_pdf)
