import numpy as np  
import ee 
import geemap 

def get_coordinates_and_utm():
    import ee
    from pyproj import Transformer
    try:
        print("Bitte geben Sie die Koordinaten als Liste ein.")
        print("Beispiel: [[lon1, lat1], [lon2, lat2], ..., [lon1, lat1]] (WGS84, also in Grad)")
        coords_input = input("Koordinaten: ")
        coords = eval(coords_input)

        polygon = ee.Geometry.Polygon(coords)
        print(f"✅ Koordinaten-Polygon in WGS84 erstellt:\n{polygon.getInfo()}")

        transformer = Transformer.from_crs("EPSG:4326", "EPSG:32632", always_xy=True)
        utm_coords = [transformer.transform(lon, lat) for lon, lat in coords]
        print("🗺️ Entsprechende UTM-Koordinaten (Zone 32N):")
        for pt in utm_coords:
            print(f"    {pt}")

        return polygon
    except Exception as e:
        print(f"❌ Fehler bei der Eingabe: {e}")
        return None
    
def get_coordinates():
    import ee
    try:
        print("Bitte geben Sie die Koordinaten als Liste ein. Beispiel: [[lon1, lat1], [lon2, lat2],         ...]")
        coords_input = input("Koordinaten: ")
        coords = eval(coords_input) 
        polygon = ee.Geometry.Polygon(coords)
        print(f"Koordinaten-Polygon erstellt: {polygon.getInfo()}")
        return polygon
    except Exception as e:
        print(f"Fehler bei der Eingabe der Koordinaten: {e}")
        return None

def calculate_mndwi_sen(image):
    """
    Computes the Modified Normalized Difference Water Index (MNDWI).

    The MNDWI is used to highlight water bodies while suppressing built-up land features.

    Formula:
        MNDWI = (Green - SWIR1) / (Green + SWIR1)

    Sentinel-2 Bands:
        - Green = B3 (560 nm)
        - SWIR1 = B11 (1610 nm)

    Parameters:
        image (ee.Image): Sentinel-2 image

    Returns:
        ee.Image: MNDWI image with values ranging from -1 to 1
    """
    import numpy as np  
    import ee 
    import geemap 
    if not isinstance(image, ee.Image):
        raise TypeError("Error: The input must be an ee.Image.")
    green = image.select('B3')
    swir = image.select('B11')
    mndwi = green.subtract(swir).divide(green.add(swir)).rename('MNDWI')
    return mndwi

def calculate_mndwi_land(image):
    """
    Computes the Modified Normalized Difference Water Index (MNDWI).

    The MNDWI is used to highlight water bodies while suppressing built-up land features.

    Formula:
        MNDWI = (Green - SWIR1) / (Green + SWIR1)

    Sentinel-2 Bands:
        - Green = B3 (560 nm)
        - SWIR1 = B11 (1610 nm)

    Parameters:
        image (ee.Image): Sentinel-2 image

    Returns:
        ee.Image: MNDWI image with values ranging from -1 to 1
    """
    import numpy as np  
    import ee 
    import geemap 
    if not isinstance(image, ee.Image):
        raise TypeError("Error: The input must be an ee.Image.")
    green = image.select('SR_B3')
    swir = image.select('SR_B6')
    mndwi = green.subtract(swir).divide(green.add(swir)).rename('MNDWI')
    return mndwi


def sunglint_removal_per_band_sen(masked_image, sample_polygon, region, bands):
    """
    Perform sunglint removal by performing a regression between SWIR1 (B11) and each band (B1 to B12)
    using the sample area, and then applying the correction globally to the MNDWI water-masked areas.

    Parameters:
        masked_image (ee.Image): The Sentinel-2 image with MNDWI-applied water mask
        sample_polygon (ee.Geometry): The user-selected deep-water sample area
        region (ee.Geometry): The entire region of interest
        bands (list): List of bands to process

    Returns:
        dict: Dictionary containing corrected band arrays
    """
    import numpy as np  # NumPy explizit importieren
    import ee  # Falls du Google Earth Engine brauchst
    import geemap  
    from sklearn.linear_model import LinearRegression
    if sample_polygon:
        sample_polygon = sample_polygon.geometry()  # Convert Feature to Geometry
        print(f"Sample polygon successfully converted to geometry.")
    else:
        print("No polygon was drawn. Please draw a polygon on the map.")
    if not sample_polygon:
        print("No sample polygon provided. Skipping sunglint correction.")
        return None

    corrected_bands = {}

    # Clip the image to the sample polygon
    sample_image = masked_image.clip(sample_polygon)

    # Select the SWIR1 (B11) band from the sample area and resample
    swir_band_sample = sample_image.select('B11').resample('bilinear').reproject(
        crs=sample_image.select('B3').projection(), scale=10
    )

    swir_array_sample = geemap.ee_to_numpy(swir_band_sample, region=sample_polygon)

    for band in bands:
        # Select and resample the current band
        current_band_sample = sample_image.select(band).resample('bilinear').reproject(
            crs=swir_band_sample.projection(), scale=10
        )

        band_array_sample = geemap.ee_to_numpy(current_band_sample, region=sample_polygon)

        # Remove NaN values and 0 values
        mask = (~np.isnan(swir_array_sample)) & (~np.isnan(band_array_sample)) & \
               (swir_array_sample != 0) & (band_array_sample != 0)

        swir_flat_sample = swir_array_sample[mask].flatten()
        band_flat_sample = band_array_sample[mask].flatten()

        # Perform linear regression
        if len(swir_flat_sample) > 0 and len(band_flat_sample) > 0:
            reg = LinearRegression().fit(swir_flat_sample[:, np.newaxis], band_flat_sample)
            slope = reg.coef_[0]

            print(f"Band: {band}, Regression slope: {slope}")

            if slope > 0:
                # Apply correction to full image
                full_swir_band = masked_image.select('B11').resample('bilinear').reproject(
                    crs=swir_band_sample.projection(), scale=10
                )
                full_band = masked_image.select(band).resample('bilinear').reproject(
                    crs=swir_band_sample.projection(), scale=10
                )

                full_swir_array = geemap.ee_to_numpy(full_swir_band, region=region)
                full_band_array = geemap.ee_to_numpy(full_band, region=region)

                corrected_band_array = full_band_array - slope * (full_swir_array - np.min(swir_flat_sample))
                corrected_band_array[corrected_band_array <= 0] = np.nan  # Avoid negative values
                
                # Restore original NaN values
                original_mask = np.isnan(full_band_array) | (full_band_array <= 0)
                corrected_band_array[original_mask] = np.nan

                corrected_bands[band] = corrected_band_array
                print(f"Band {band}: Sunglint correction applied (slope = {slope:.4f}).")
            else:
                print(f"Band {band}: Negative slope. Using original MNDWI-masked band.")
                corrected_bands[band] = geemap.ee_to_numpy(masked_image.select(band), region=region)
        else:
            print(f"Band {band}: Not enough valid pixels. Using original MNDWI-masked band.")
            corrected_bands[band] = geemap.ee_to_numpy(masked_image.select(band), region=region)

    return corrected_bands

def sunglint_removal_per_band_land(masked_image, sample_polygon, region, bands):
    """
    Perform sunglint removal by performing a regression between SWIR1 (B11) and each band (B1 to B12)
    using the sample area, and then applying the correction globally to the MNDWI water-masked areas.

    Parameters:
        masked_image (ee.Image): The Sentinel-2 image with MNDWI-applied water mask
        sample_polygon (ee.Geometry): The user-selected deep-water sample area
        region (ee.Geometry): The entire region of interest
        bands (list): List of bands to process

    Returns:
        dict: Dictionary containing corrected band arrays
    """
    import numpy as np  # NumPy explizit importieren
    import ee  # Falls du Google Earth Engine brauchst
    import geemap  
    from sklearn.linear_model import LinearRegression
    if sample_polygon:
        sample_polygon = sample_polygon.geometry()  # Convert Feature to Geometry
        print(f"Sample polygon successfully converted to geometry.")
    else:
        print("No polygon was drawn. Please draw a polygon on the map.")
    if not sample_polygon:
        print("No sample polygon provided. Skipping sunglint correction.")
        return None

    corrected_bands = {}

    # Clip the image to the sample polygon
    sample_image = masked_image.clip(sample_polygon)

    # Select the SWIR1 (B11) band from the sample area and resample
    swir_band_sample = sample_image.select('SR_B6').resample('bilinear').reproject(
        crs=sample_image.select('SR_B3').projection(), scale=10
    )

    swir_array_sample = geemap.ee_to_numpy(swir_band_sample, region=sample_polygon)

    for band in bands:
        # Select and resample the current band
        current_band_sample = sample_image.select(band).resample('bilinear').reproject(
            crs=swir_band_sample.projection(), scale=10
        )

        band_array_sample = geemap.ee_to_numpy(current_band_sample, region=sample_polygon)

        # Remove NaN values and 0 values
        mask = (~np.isnan(swir_array_sample)) & (~np.isnan(band_array_sample)) & \
               (swir_array_sample != 0) & (band_array_sample != 0)

        swir_flat_sample = swir_array_sample[mask].flatten()
        band_flat_sample = band_array_sample[mask].flatten()

        # Perform linear regression
        if len(swir_flat_sample) > 0 and len(band_flat_sample) > 0:
            reg = LinearRegression().fit(swir_flat_sample[:, np.newaxis], band_flat_sample)
            slope = reg.coef_[0]

            print(f"Band: {band}, Regression slope: {slope}")

            if slope > 0:
                # Apply correction to full image
                full_swir_band = masked_image.select('SR_B6').resample('bilinear').reproject(
                    crs=swir_band_sample.projection(), scale=10
                )
                full_band = masked_image.select(band).resample('bilinear').reproject(
                    crs=swir_band_sample.projection(), scale=10
                )

                full_swir_array = geemap.ee_to_numpy(full_swir_band, region=region)
                full_band_array = geemap.ee_to_numpy(full_band, region=region)

                corrected_band_array = full_band_array - slope * (full_swir_array - np.min(swir_flat_sample))
                corrected_band_array[corrected_band_array <= 0] = np.nan  # Avoid negative values
                
                # Restore original NaN values
                original_mask = np.isnan(full_band_array) | (full_band_array <= 0)
                corrected_band_array[original_mask] = np.nan

                corrected_bands[band] = corrected_band_array
                print(f"Band {band}: Sunglint correction applied (slope = {slope:.4f}).")
            else:
                print(f"Band {band}: Negative slope. Using original MNDWI-masked band.")
                corrected_bands[band] = geemap.ee_to_numpy(masked_image.select(band), region=region)
        else:
            print(f"Band {band}: Not enough valid pixels. Using original MNDWI-masked band.")
            corrected_bands[band] = geemap.ee_to_numpy(masked_image.select(band), region=region)

    return corrected_bands

def calculate_subsurface_angle(theta_deg, nW=1.33):
        """Applies Snell's Law to compute subsurface angles (returns degrees)."""
        import ee
        import numpy as np
        theta_rad = np.radians(theta_deg)  # Convert to radians
        theta_sub_rad = np.arcsin(np.sin(theta_rad) / nW)  # Apply Snell's Law
        return np.degrees(theta_sub_rad)  # Convert back to degrees
    
def extract_and_convert_metadata(image, bands=["B1","B2","B3", "B4","B5", "B6","B7", "B11"], nW=1.33):
    """
    Extracts metadata for viewing angles, solar zenith angles, and converts them to subsurface angles.
    1. Converts angles using Snell’s Law (keeps degrees).
    2. Converts all angles into radians after Snell's Law.

    Parameters:
        image (ee.Image): Sentinel-2 image
        bands (list): List of bands to extract viewing angles for
        nW (float): Water refractive index (default = 1.33 for freshwater)

    Returns:
        dict: Metadata dictionary with converted angles
    """
    import ee
    import numpy as np
    def calculate_subsurface_angle(theta_deg, nW=1.33):
        """Applies Snell's Law to compute subsurface angles (returns degrees)."""
        theta_rad = np.radians(theta_deg)  # Convert to radians
        theta_sub_rad = np.arcsin(np.sin(theta_rad) / nW)  # Apply Snell's Law
        return np.degrees(theta_sub_rad)  # Convert back to degrees

    metadata = {
        'DateTime': ee.Date(image.get('system:time_start')).format("YYYY-MM-dd HH:mm:ss").getInfo(),
        'Solar Zenith Angle (θs) (deg)': image.get('MEAN_SOLAR_ZENITH_ANGLE').getInfo()
    }

    # Compute subsurface solar zenith angle in degrees
    metadata['Subsurface Solar Zenith Angle (θs) (deg)'] = calculate_subsurface_angle(metadata['Solar Zenith Angle (θs) (deg)'], nW)
    
    # Convert both angles to radians
    metadata['Solar Zenith Angle (θs) (rad)'] = np.radians(metadata['Solar Zenith Angle (θs) (deg)'])
    metadata['Subsurface Solar Zenith Angle (θs) (rad)'] = np.radians(metadata['Subsurface Solar Zenith Angle (θs) (deg)'])

    for band in bands:
        key = f'Mean Viewing Angle {band} (θv) (deg)'
        property_name = f'MEAN_INCIDENCE_ZENITH_ANGLE_{band}'
        theta_v_deg = image.get(property_name).getInfo()

        metadata[key] = theta_v_deg  # Store original value
        metadata[f'Subsurface Viewing Angle {band} (θv) (deg)'] = calculate_subsurface_angle(theta_v_deg, nW)

        # Convert both to radians
        metadata[f'Mean Viewing Angle {band} (θv) (rad)'] = np.radians(metadata[key])
        metadata[f'Subsurface Viewing Angle {band} (θv) (rad)'] = np.radians(metadata[f'Subsurface Viewing Angle {band} (θv) (deg)'])

    return metadata


def calc_w(a, bb):
    """Computes w as the ratio of backscatter (bb) to total attenuation (a + bb)."""
    return bb / (a + bb)

def calculate_Kd(ko,aW, bW, theta_s):
    import numpy as np
    """
    Computes Kd (diffuse attenuation coefficient).
    
    Formula:
        Kd = coef1 * (aW + bW) / cos(theta_s)
    
    Parameters:
        aW (float): Absorption coefficient
        bW (float): Backscatter coefficient
        theta_s (float): Subsurface solar zenith angle in radians
    
    Returns:
        float: Kd value
    """
    return ko * (aW + bW) / np.cos(theta_s)

def calculate_Rrs_infinity(prs1, prs2,prs3,prs4,prs5,prs6,prs7,theta_v, theta_s, w):
    """
    Computes Rrs_infinity using an empirical polynomial model.
    
    Formula:
        Rrs_infinity = coef1 * (1 + coef2 * w + coef3 * w^2 + ...) * 
                       (1 + coefX * (1 / cos(theta_s))) * (1 + coefY * (1 / cos(theta_v))) * w
    
    Parameters:
        theta_v (float): Subsurface viewing angle in radians
        theta_s (float): Subsurface solar zenith angle in radians
        w (float): Computed value from calc_w()
    
    Returns:
        float: Rrs_infinity value
    """
    import numpy as np
    Rrs_infinity = prs1 * (1 + prs2 * w + prs3 * w**2 + prs4 * w**3) * \
                   (1 + (prs5 * (1 / np.cos(theta_s)))) * \
                   (1+ prs6*0) * \
                   (1 + (prs7 * (1 / np.cos(theta_v))) )*w
    return Rrs_infinity

def calculate_f_arrow(image_idx, results_w, rrs_infinity_results):
    """
    Berechnet f_arrow für ein gegebenes Bild basierend auf Band B7.

    Parameter:
        image_idx (int): Die Bildnummer (Index beginnt bei 1)
        results_w (list): Liste der berechneten w-Werte für alle Bilder
        rrs_infinity_results (list): Liste der berechneten Rrs_infinity-Werte für alle Bilder

    Rückgabe:
        float oder None: Berechneter f_arrow-Wert oder None, falls Daten fehlen
    """
    import numpy as np
    w_value = next((r['w'] for r in results_w if r['Image'] == image_idx and r['Band'] == 'B7'), None)
    
    rrs_inf = next((r['Rrs_infinity'] for r in rrs_infinity_results if r['Image'] == image_idx and r['Band'] == 'B7'), None)

    if w_value is None or rrs_inf is None:
        print(f"Skipping Image {image_idx} due to missing values for Band B7.")
        return None

    f_arrow = rrs_inf / w_value

    return f_arrow

def calculate_R_u(R_rs, R_B, K_d, theta_v, z_B, f_u, ars1, ars2):
    """
    Berechnet das aufwärtsgerichtete Reflektanzsignal `R_u`.

    Parameter:
        R_rs (array): Gemessene Wasserreflektanz
        R_B (float): Bodenreflektanz
        K_d (float): Diffuser Lichtausbreitungskoeffizient
        theta_v (float): Subsurface Viewing Angle in Radiant
        z_B (array): Wasser-Tiefenwerte (aus Mosaik)
        f_u (float): `f_arrow`-Wert für das Bild
        ars1 (float): Empirischer Koeffizient für den oberen Teil der Formel
        ars2 (float): Empirischer Koeffizient für den unteren Teil der Formel

    Rückgabe:
        np.array: Berechnetes `R_u`
    """
    import numpy as np
    
    exp_term = np.exp(-K_d * (1 + 1 / np.cos(theta_v)) * z_B)

    R_u = (R_rs - ars2 * (R_B / np.pi) * exp_term) / (f_u * (1 - ars1 * exp_term))

    return R_u

def calculate_Cx(R_u, aw, bbw, bb_x):
    """
    Berechnet die Konzentration C_X basierend auf R_u.

    Parameter:
        R_u (array): Aufwärtsgerichtete Reflektanz R_u
        aw (float): Absorptionskoeffizient des Wassers für das Band
        bbw (float): Rückstreukoeffizient für das Band

    Rückgabe:
        np.array: Berechneter C_X-Wert
    """
    import numpy as np
    numerator = R_u * (aw + bbw) - bbw
    denominator = bb_x * (1 - R_u)  

    C_X = numerator / denominator

    return C_X

def calculate_Ku(a, bb, wB, theta_s, k1, k2):
    import numpy as np
    """
    Berechnet den Aufwärts-Diffusionsdämpfungskoeffizienten Ku.
    
    Parameter:
        a (float): Absorptionskoeffizient
        bb (float): Rückstreukoeffizient
        wB (float): Single-Scattering-Albedo
        theta_s (float): Solar Zenith Angle (Subsurface)
        k1 (float): Empirischer Parameter k1
        k2 (float): Empirischer Parameter k2

    Rückgabe:
        float: Ku-Wert
    """
    return (a + bb) * (1 + wB) ** k1 * (1 + k2 * (1 / np.cos(theta_s)))

def optimize_a(R_rs_measured, C_X, empirical_coefficients, parameters, theta_s, theta_v, R_B, z_B, band, k1u, k2u, k1b, k2b, ars1, ars2, ko, a_0):
    """
    Iteriert über `a` für jeden Pixel, um den Fehler zwischen `R_rs_measured` und `R_rs_simulated` zu minimieren.
    NaN-Werte bleiben erhalten und beeinflussen die Berechnung nicht.
    """
    a_0 = 2
    imax = 5000
    threshold = 0.01
    import numpy as np
    R_rs_measured = np.squeeze(R_rs_measured)  # Entfernt unnötige Dimensionen
    C_X = np.squeeze(C_X)
    z_B = np.squeeze(z_B)
    R_B = np.squeeze(R_B)

    R_rs_measured = np.array(R_rs_measured, dtype=np.float32)
    C_X = np.array(C_X, dtype=np.float32)
    z_B = np.array(z_B, dtype=np.float32)
    
    # 🔹 Maske für gültige Werte (keine NaNs)
    valid_mask = ~np.isnan(R_rs_measured) & ~np.isnan(C_X) & ~np.isnan(z_B)

    # Falls keine gültigen Werte vorhanden sind → direkt `NaN` zurückgeben
    if np.sum(valid_mask) == 0:
        print("❌ Keine gültigen Werte vorhanden, überspringe Optimierung.")
        return np.full_like(R_rs_measured, np.nan)

    # 🔹 `a` initialisieren, aber nur für gültige Pixel
    a = np.full_like(R_rs_measured, np.nan)
    a[valid_mask] = a_0  # Nur wo `valid_mask` True ist, wird `a_0` gesetzt

    # 🔹 Mindestwert für `a` setzen
    aw = empirical_coefficients[band]["aw"]
    
    for iteration in range(imax):
        # Berechnung nur für gültige Pixel
        bbw = empirical_coefficients[band]["bbw"]
        bb_x = empirical_coefficients[band]["b*b_X"]
        bb = np.full_like(a, np.nan)
        bb[valid_mask] = bbw + (bb_x * C_X[valid_mask])

        # Berechnung von `wB`
        wB = np.full_like(a, np.nan)
        wB[valid_mask] = bb[valid_mask] / (a[valid_mask] + bb[valid_mask])

        # Berechnung von `Kd`
        Kd = np.full_like(a, np.nan)
        Kd[valid_mask] = ko * (a[valid_mask] + bb[valid_mask]) / np.cos(theta_s)

        # Berechnung von `Rrs_∞`
        Rrs_inf = np.full_like(a, np.nan)
        Rrs_inf[valid_mask] = ars1 * (1 + ars2 * wB[valid_mask] + ars2 * wB[valid_mask]**2) * \
                              (1 + (ars2 * (1 / np.cos(theta_s)))) * (1 + (ars2 * (1 / np.cos(theta_v)))) * wB[valid_mask]

        # Berechnung von `KuW` und `KuB`
        Ku_W = np.full_like(a, np.nan)
        Ku_B = np.full_like(a, np.nan)
        Ku_W[valid_mask] = (a[valid_mask] + bb[valid_mask]) * (1 + wB[valid_mask]) ** k1u * \
                           (1 + k2u * (1 / np.cos(theta_s))) / np.cos(theta_v)
        Ku_B[valid_mask] = (a[valid_mask] + bb[valid_mask]) * (1 + wB[valid_mask]) ** k1b * \
                           (1 + k2b * (1 / np.cos(theta_s))) / np.cos(theta_v)

        # Berechnung von `R_rs_simulated`
        exp_KuW = np.exp(- (Kd+Ku_W) * z_B)
        exp_KuB = np.exp(- (Kd+Ku_B) * z_B)

        R_rs_simulated = np.full_like(a, np.nan)
        R_rs_simulated[valid_mask] = Rrs_inf[valid_mask] * (1 - ars1 * exp_KuW[valid_mask]) + \
                                     ars2 * (R_B / np.pi) * exp_KuB[valid_mask]

        # Fehlerberechnung (nur gültige Werte)
        error = np.full_like(a, np.nan)
        error[valid_mask] = np.abs(R_rs_measured[valid_mask] - R_rs_simulated[valid_mask])

        # Falls Fehler zu groß, iteriere weiter
        if np.nanmean(error) < threshold:
            print(f"✅ Konvergenz nach {iteration} Iterationen erreicht.")
            return a

        # 🔹 `a` aktualisieren, aber nur für gültige Werte
        a[valid_mask] = np.maximum(a[valid_mask] * 0.99 + R_rs_measured[valid_mask] * 0.01, aw)

    print(f"⚠️ Maximale Iterationen erreicht. Letzter `a`-Wert: Mittelwert = {np.nanmean(a):.4f}")
    return a

def calculate_cp_ay(a_measured, a_W_values, a_P_star_values, lambda_values, s_Y=0.011, lambda_0=400):
    """
    Berechnet die unbekannten Parameter C_P und a_Y(lambda_0) aus dem gegebenen Gleichungssystem.

    Parameter:
    a_measured : array
        Gemessene Gesamtabsorption a(lambda) für verschiedene Wellenlängen.
    a_W_values : array
        Absorptionskoeffizienten des Wassers für die entsprechenden Wellenlängen.
    a_P_star_values : array
        Spezifische Absorptionskoeffizienten der Partikel für die entsprechenden Wellenlängen.
    lambda_values : array
        Wellenlängen (nm), für die die Berechnung durchgeführt wird.
    s_Y : float, optional
        Exponentielle Abfallrate, standardmäßig 0.011 nm^-1.
    lambda_0 : int, optional
        Referenzwellenlänge in nm, standardmäßig 400 nm.

    Rückgabe:
    (C_P, a_Y_lambda0) : tuple
        Optimierte Werte für C_P und a_Y(lambda_0).
    """
    import numpy as np
    from scipy.optimize import least_squares

    def residuals(params):
        C_P, a_Y_lambda0 = params
        return a_measured - (a_W_values + a_P_star_values * C_P + a_Y_lambda0 * np.exp(-s_Y * (lambda_values - lambda_0)))

    # Startwerte für die Optimierung
    initial_guess = [10, 0.5]  # Startwerte für C_P und a_Y(lambda0)

    # Least Squares Optimierung zur Lösung des Gleichungssystems
    result = least_squares(residuals, initial_guess)

    return result.x

def load_existing_depth(filepath, expected_extent, expected_crs, expected_resolution=(10, 10)):
    import os
    import rasterio
    from rasterio.enums import Resampling
    import numpy as np
    
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return None, None

    with rasterio.open(filepath) as src:
        data = src.read(1)  # assumes single-band depth data
        crs = src.crs
        transform = src.transform
        bounds = src.bounds
        resolution = (src.res[0], src.res[1])

        if crs != expected_crs:
            print(f"Reprojecting {filepath} to match expected CRS...")
            dst_transform, width, height = calculate_default_transform(
                src.crs, expected_crs, src.width, src.height, *src.bounds)

            kwargs = src.meta.copy()
            kwargs.update({
                'crs': expected_crs,
                'transform': dst_transform,
                'width': width,
                'height': height
            })

            reprojected = np.empty((height, width), dtype=np.float32)

            reproject(
                source=rasterio.band(src, 1),
                destination=reprojected,
                src_transform=transform,
                src_crs=src.crs,
                dst_transform=dst_transform,
                dst_crs=expected_crs,
                resampling=Resampling.bilinear
            )

            transform = dst_transform
            bounds = rasterio.transform.array_bounds(height, width, transform)
            resolution = (transform.a, -transform.e)  # pixel size
            data = reprojected

        else:
            print(f"Using original CRS ({expected_crs})")
            data = src.read(1)

        src_extent = [bounds.left, bounds.right, bounds.bottom, bounds.top]
        lon_overlap = not (src_extent[1] <= expected_extent[0] or src_extent[0] >= expected_extent[1])
        lat_overlap = not (src_extent[3] <= expected_extent[2] or src_extent[2] >= expected_extent[3])
        
        if not (lon_overlap and lat_overlap):
            print(f"Skipping {filepath}: No overlap with expected region.")
            return None, None

        if resolution != expected_resolution:
            print(f"Interpolating {filepath} to {expected_resolution} m resolution...")
            scale_x = resolution[0] / expected_resolution[0]
            scale_y = resolution[1] / expected_resolution[1]
            data = src.read(
                1,
                out_shape=(
                    int(src.height * scale_y),
                    int(src.width * scale_x)
                ),
                resampling=Resampling.bilinear
            )

        return data, src_extent
    
def calculate_subsurface_angle(theta_deg, nW=1.33):
        import numpy as np
        theta_rad = np.radians(theta_deg)
        theta_sub_rad = np.arcsin(np.sin(theta_rad) / nW)
        return np.degrees(theta_sub_rad)
    
def extract_metadata_solar_only(image, bands=["SR_B1", "SR_B2", "SR_B3", "SR_B4", "SR_B5", "SR_B6"], nW=1.33):
        import numpy as np
        import ee

        try:
            timestamp = ee.Date(image.get("system:time_start")).format("YYYY-MM-dd HH:mm:ss").getInfo()
            solar_elev = image.get("SUN_ELEVATION").getInfo()
            solar_zenith_deg = float(90 - solar_elev)
            subsurface_sza_deg = calculate_subsurface_angle(solar_zenith_deg, nW)

            metadata = {
                "Timestamp": timestamp,
                "Solar Angle (deg)": solar_zenith_deg,
                "Subsurface Solar Angle (deg)": subsurface_sza_deg,
                "Solar Angle (rad)": np.radians(solar_zenith_deg),
                "Subsurface Solar Angle (rad)": np.radians(subsurface_sza_deg),
                }
                
            for band in bands:
        
                metadata[f"Mean Viewing Angle {band} (θv) (deg)"] = 0.0
                metadata[f"Subsurface Viewing Angle {band} (θv) (deg)"] = 0.0
                metadata[f"Mean Viewing Angle {band} (θv) (rad)"] = 0.0
                metadata[f"Subsurface Viewing Angle {band} (θv) (rad)"] = 0.0

            return metadata

        except Exception as e:
            return {
                "Timestamp": "Fehler",
                "Fehlermeldung": str(e)
            }
        
def match_reference_tiff_to_band_depth(tiff_path, extent_b7, shape_b7, b7_crs="EPSG:4326", tiff_crs="EPSG:32632", label=""):
    from pyproj import Transformer
    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.warp import Resampling
    import matplotlib as plt
    if len(shape_b7) == 3:
        height, width = shape_b7[0:2]
    elif len(shape_b7) == 2:
        height, width = shape_b7
    else:
        raise ValueError("❌ Ungültiges Format für shape_b7")

    transformer = Transformer.from_crs(b7_crs, tiff_crs, always_xy=True)
    x_min, y_min = transformer.transform(extent_b7[0], extent_b7[2])
    x_max, y_max = transformer.transform(extent_b7[1], extent_b7[3])

    with rasterio.open(tiff_path) as src:
        try:
            window = from_bounds(x_min, y_min, x_max, y_max, transform=src.transform)
            raw_data = src.read(
                1,
                window=window,
                out_shape=(height, width),
                resampling=Resampling.bilinear
            )

            fill_value = 0
            if src.nodata is not None:
                raw_data = np.where(raw_data == src.nodata, fill_value, raw_data)

            src_data = np.nan_to_num(raw_data, nan=fill_value)
            src_data[src_data == 0] = np.nan

        except Exception as e:
            print(f"❌ Fehler beim Lesen & Interpolieren ({label}): {e}")
            return None

    plt.figure(figsize=(10, 6))
    plt.imshow(src_data, cmap=cmocean.cm.deep, extent=extent_b7, origin='upper')
    plt.title(f"{label} (zugeschnitten & angepasst an B7)")
    plt.xlabel("Longitude [°]")
    plt.ylabel("Latitude [°]")
    plt.colorbar(label=f"Depth [m]")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return src_data
def calculate_depth_Albert(ars1, ars2, r_rs_measured, r_rs_inf, R_B, Kd, Ov):
    """
    Berechnung der Wassertiefe basierend auf Albert et al.
    
    Parameter:
    r_rs_measured : float oder np.array
        Gemessene Reflektanz unter Wasser
    r_rs_inf : float
        Theoretische Reflektanz für unendliche Wassertiefe
    R_B : float
        Bodenreflektanz
    Kd : float
        Diffuser Lichtausbreitungskoeffizient
    Ov : float
        Unterwasser-Betrachtungswinkel
    
    Rückgabe:
    h : np.array
        Berechnete Wassertiefe in Metern
    """
    import numpy as np
    import matplotlib.pyplot as plt
    
    
    ln_argument= ((ars1 * r_rs_inf- ars2*(R_B/np.pi)/(r_rs_inf-r_rs_measured)))
    h =  (1 / (Kd * (1 + (1 / np.cos(Ov))))) * np.log(ln_argument)
    return h

def plot_comparison(name_x, data_x, name_y, data_y, min_depth=None, max_depth=None):
    valid_mask = ~np.isnan(data_x) & ~np.isnan(data_y)

    if max_depth is not None:
        valid_mask &= data_y <= max_depth
    if min_depth is not None:
        valid_mask &= data_y >= min_depth

    x = data_x[valid_mask].flatten()
    y = data_y[valid_mask].flatten()

    if len(x) == 0:
        print(f"⚠️ Keine gültigen Daten für {name_x} vs. {name_y} im Bereich {min_depth}–{max_depth} m.")
        return
    
def evaluate_bins(name_x, data_x, name_y, data_y, bin_size=0.5):
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
    import pandas as pd
    valid_mask = ~np.isnan(data_x) & ~np.isnan(data_y)
    x = data_x[valid_mask].flatten()
    y = data_y[valid_mask].flatten()

    max_depth = np.ceil(np.max(y))
    bin_edges = np.arange(0, max_depth + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    stats = {
        "Tiefenbereich [m]": [],
        "N": [],
        "RMSE [m]": [],
        "MAE [m]": [],
        "Bias [m]": [],
        "R²": [],
    }

    for i in range(len(bin_edges) - 1):
        mask = (y >= bin_edges[i]) & (y < bin_edges[i + 1])
        if np.sum(mask) > 1:
            stats["Tiefenbereich [m]"].append(f"{bin_edges[i]:.1f}–{bin_edges[i+1]:.1f}")
            stats["N"].append(np.sum(mask))
            stats["RMSE [m]"].append(mean_squared_error(y[mask], x[mask], squared=False))
            stats["MAE [m]"].append(mean_absolute_error(y[mask], x[mask]))
            stats["Bias [m]"].append(np.mean(x[mask] - y[mask]))
            stats["R²"].append(r2_score(y[mask], x[mask]))
        else:
            stats["Tiefenbereich [m]"].append(f"{bin_edges[i]:.1f}–{bin_edges[i+1]:.1f}")
            stats["N"].append(0)
            stats["RMSE [m]"].append(np.nan)
            stats["MAE [m]"].append(np.nan)
            stats["Bias [m]"].append(np.nan)
            stats["R²"].append(np.nan)

    df = pd.DataFrame(stats)
    print(f"\n📊 Statistik-Tabelle: {name_x} vs. {name_y}")
    print(df.to_string(index=False))
    return df

def match_reference_tiff_to_band_cx(tiff_path, extent_b7, shape_b7, b7_crs="EPSG:4326", tiff_crs="EPSG:32632", label=""):
    from pyproj import Transformer
    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.warp import Resampling
    import matplotlib.pyplot as plt
    
    if len(shape_b7) == 3:
        height, width = shape_b7[0:2]
    elif len(shape_b7) == 2:
        height, width = shape_b7
    else:
        raise ValueError("❌ Ungültiges Format für shape_b7")

    transformer = Transformer.from_crs(b7_crs, tiff_crs, always_xy=True)
    x_min, y_min = transformer.transform(extent_b7[0], extent_b7[2])
    x_max, y_max = transformer.transform(extent_b7[1], extent_b7[3])

    with rasterio.open(tiff_path) as src:
        try:
            window = from_bounds(x_min, y_min, x_max, y_max, transform=src.transform)
            raw_data = src.read(
                1,
                window=window,
                out_shape=(height, width),
                resampling=Resampling.bilinear
            )

            fill_value = 0
            if src.nodata is not None:
                raw_data = np.where(raw_data == src.nodata, fill_value, raw_data)

            src_data = np.nan_to_num(raw_data, nan=fill_value)
            src_data[src_data == 0] = np.nan

        except Exception as e:
            print(f"❌ Fehler beim Lesen & Interpolieren ({label}): {e}")
            return None

    
    plt.figure(figsize=(10, 6))
    plt.imshow(src_data, cmap=cmocean.cm.matter, extent=extent_b7, origin='upper')
    plt.title(f"{label} (zugeschnitten & angepasst an B7)")
    plt.xlabel("Longitude [°]")
    plt.ylabel("Latitude [°]")
    plt.colorbar(label=f"Suspended Matter Concentration [mg/l]")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return src_data

def calculate_f_arrow_land(image_idx, results_w, rrs_infinity_results):
                """
                Berechnet f_arrow für ein gegebenes Bild basierend auf Band B7.

                Parameter:
                    image_idx (int): Die Bildnummer (Index beginnt bei 1)
                    results_w (list): Liste der berechneten w-Werte für alle Bilder
                    rrs_infinity_results (list): Liste der berechneten Rrs_infinity-Werte für alle Bilder

                Rückgabe:
                    float oder None: Berechneter f_arrow-Wert oder None, falls Daten fehlen
                """
                import numpy as np
          
                w_value = next((r['w'] for r in results_w if r['Image'] == image_idx and r['Band'] == 'SR_B5'), None)

                rrs_inf = next((r['Rrs_infinity'] for r in rrs_infinity_results if r['Image'] == image_idx and r['Band'] == 'SR_B5'), None)

                if w_value is None or rrs_inf is None:
                    print(f"Skipping Image {image_idx} due to missing values for Band B7.")
                    return None
                f_arrow = rrs_inf / w_value

                return f_arrow
            
def plot_comparison(name_x, data_x, name_y, data_y, max_depth=None):
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
    import matplotlib as plt
    valid_mask = ~np.isnan(data_x) & ~np.isnan(data_y)
    if max_depth is not None:
        valid_mask &= data_y <= max_depth

    x = data_x[valid_mask].flatten()
    y = data_y[valid_mask].flatten()

    if len(x) == 0:
        print(f"⚠️ Keine gültigen Daten für {name_x} vs. {name_y} bis {max_depth} m.")
        return
    r2 = r2_score(y, x)
    rmse = mean_squared_error(y, x, squared=False)
    mae = mean_absolute_error(y, x)
    bias = np.mean(x - y)
    
    coeffs = np.polyfit(x, y, 1)
    poly_eqn = np.poly1d(coeffs)
    x_line = np.linspace(min(x.min(), y.min()), max(x.max(), y.max()), 100)
    y_line = poly_eqn(x_line)

    distance = np.abs(x - y)

    fig, ax = plt.subplots(figsize=(6.5, 5)) 
    sc = ax.scatter(x, y, c=distance, cmap='hot', s=10, alpha=0.6)
    ax.plot(x_line, x_line, color='black', linestyle='--', label='1:1 Linie')
    ax.plot(x_line, y_line, color='red', label='Regressionslinie')

    ax.set_xlabel(f"{name_x} depth [m]")
    ax.set_ylabel(f"{name_y} depth [m]")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right")

    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Abstand zur 1:1-Linie [m]")

    stats_text = (
        f"$N$ = {len(x)}\n"
        f"$R^2$ = {r2:.2f}\n"
        f"RMSE = {rmse:.2f} m\n"
        f"MAE = {mae:.2f} m\n"
        f"Bias = {bias:.2f} m"
    )
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    plt.tight_layout()
    plt.show()
    
def match_reference_tiff_to_band_cp(tiff_path, extent_b7, shape_b7, b7_crs="EPSG:4326",  tiff_crs="EPSG:32632", label=""):
    from pyproj import Transformer
    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.warp import Resampling
    import matplotlib.pyplot as plt
    if len(shape_b7) == 3:
        height, width = shape_b7[0:2]
    elif len(shape_b7) == 2:
        height, width = shape_b7
    else:
        raise ValueError("❌ Ungültiges Format für shape_b7")

    transformer = Transformer.from_crs(b7_crs, tiff_crs, always_xy=True)
    x_min, y_min = transformer.transform(extent_b7[0], extent_b7[2])
    x_max, y_max = transformer.transform(extent_b7[1], extent_b7[3])

    with rasterio.open(tiff_path) as src:
        try:
            window = from_bounds(x_min, y_min, x_max, y_max, transform=src.transform)
            raw_data = src.read(
                1,
                window=window,
                out_shape=(height, width),
                resampling=Resampling.bilinear
            )

            fill_value = 0
            if src.nodata is not None:
                raw_data = np.where(raw_data == src.nodata, fill_value, raw_data)

            src_data = np.nan_to_num(raw_data, nan=fill_value)
            src_data[src_data == 0] = np.nan

        except Exception as e:
            print(f"❌ Fehler beim Lesen & Interpolieren ({label}): {e}")
            return None

    # Darstellung
    plt.figure(figsize=(10, 6))
    plt.imshow(src_data, cmap=cmocean.cm.algae, extent=extent_b7, origin='upper')
    plt.title(f"{label} (zugeschnitten & angepasst an B7)")
    plt.xlabel("Longitude [°]")
    plt.ylabel("Latitude [°]")
    plt.colorbar(label=f"Phytoplancton Concentration [µg/l]")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return src_data

def match_reference_tiff_to_band_ay(tiff_path, extent_b7, shape_b7, b7_crs="EPSG:4326", tiff_crs="EPSG:32632", label=""):
    from pyproj import Transformer
    import rasterio
    from rasterio.windows import from_bounds
    from rasterio.warp import Resampling
    import matplotlib.pyplot as plt
    if len(shape_b7) == 3:
        height, width = shape_b7[0:2]
    elif len(shape_b7) == 2:
        height, width = shape_b7
    else:
        raise ValueError("❌ Ungültiges Format für shape_b7")

    transformer = Transformer.from_crs(b7_crs, tiff_crs, always_xy=True)
    x_min, y_min = transformer.transform(extent_b7[0], extent_b7[2])
    x_max, y_max = transformer.transform(extent_b7[1], extent_b7[3])

    with rasterio.open(tiff_path) as src:
        try:
            window = from_bounds(x_min, y_min, x_max, y_max, transform=src.transform)
            raw_data = src.read(
                1,
                window=window,
                out_shape=(height, width),
                resampling=Resampling.bilinear
            )

            fill_value = 0
            if src.nodata is not None:
                raw_data = np.where(raw_data == src.nodata, fill_value, raw_data)

            src_data = np.nan_to_num(raw_data, nan=fill_value)
            src_data[src_data == 0] = np.nan

        except Exception as e:
            print(f"❌ Fehler beim Lesen & Interpolieren ({label}): {e}")
            return None

    # Darstellung
    plt.figure(figsize=(10, 6))
    plt.imshow(src_data, cmap=cmocean.cm.turbid, extent=extent_b7, origin='upper')
    plt.title(f"{label} (zugeschnitten & angepasst an B7)")
    plt.xlabel("Longitude [°]")
    plt.ylabel("Latitude [°]")
    plt.colorbar(label=f"( a_Y(lambda_0) [m^-1]")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

    return src_data

def simulate_rrs_pixel(a, theta_s, theta_v, Rrs_inf, ars1, ars2, k0,
                       krs1_W, krs2_W, krs1_B, krs2_B, R_B, z_B, bb):
    import numpy as np
    wb = bb / (a + bb)
    
    Kd = k0 * (a + bb) / np.cos(theta_s)
    Ku_W = (a + bb) * (1 + wb)**krs1_W * (1 + krs2_W / np.cos(theta_s)) / np.cos(theta_v)
    Ku_B = (a + bb) * (1 + wb)**krs1_B * (1 + k2b / np.cos(theta_s)) / np.cos(theta_v)

    exp1 = np.exp(-(Kd + Ku_W) * z_B)
    exp2 = np.exp(-(Kd + Ku_B) * z_B)

    Rrs_sim = Rrs_inf * (1 - ars1 * exp1) + ars2 * (R_B / np.pi) * exp2
    return Rrs_sim

def residuals_pixelwise(a, Rrs_measured, theta_s, theta_v, Rrs_inf,
                        ars1, ars2, k0, krs1_W, krs2_W, krs1_B, krs2_B,
                        R_B, z_B, bb):
    Rrs_model = simulate_rrs_pixel(a[0], theta_s, theta_v, Rrs_inf, ars1, ars2,
                                   k0, krs1_W, krs2_W, krs1_B, krs2_B, R_B, z_B, bb)
    return Rrs_measured - Rrs_model

def optimize_a_pixelwise(
    R_rs_measured, theta_s, theta_v, Rrs_inf,
    ars1, ars2, k0, k1u, k2u, k1b, k2b,
    R_B, z_B, bb, aw=1.0, a_init=0.75, threshold=0.001, max_iter=2000
):
    import numpy as np

    a_result = np.full_like(R_rs_measured, np.nan, dtype=np.float64)
    valid_mask = ~np.isnan(R_rs_measured) & ~np.isnan(z_B) & ~np.isnan(bb)
    rows, cols = np.where(valid_mask)

    for i, j in zip(rows, cols):
        a = a_init
        for _ in range(max_iter):
            wb = bb[i, j] / (a + bb[i, j])
            
            Kd = k0 * (a + bb[i, j]) / np.cos(theta_s)
            Ku_W = ((a + bb[i, j]) * (1 + wb)**k1u * (1 + k2u / np.cos(theta_s)))/ np.cos(theta_v)
            Ku_B = ((a + bb[i, j]) * (1 + wb)**k1b * (1 + k2b / np.cos(theta_s)))/ np.cos(theta_v)

            exp1 = np.exp(-(Kd + Ku_W) * z_B[i, j])
            exp2 = np.exp(-(Kd + Ku_B) * z_B[i, j])

            Rrs_sim = Rrs_inf[i, j] * (1 - ars1 * exp1) + ars2 * (R_B / np.pi) * exp2
            error = abs(R_rs_measured[i, j] - Rrs_sim)

            if error < threshold:
                break

            a += (R_rs_measured[i, j] - Rrs_sim) * 0.005
            a = np.clip(a, aw, 5.0)

        a_result[i, j] = a

    return a_result

def residuals(params, a_measured, a_W, a_P_star, lambda_val):
        C_P, a_Y_lambda0 = params
        return a_measured - (a_W + a_P_star * C_P + a_Y_lambda0 * np.exp(-0.011 * (lambda_val - LAMBDA_0)))
    
def residuals_opt(params, Rrs_pixel, theta_s, theta_v, Rrs_inf, ars1, ars2, ko, aw, aP_star, wavelength, bbw, bb_x):
    return Rrs_pixel - simulate_rrs_pixel(params, theta_s, theta_v, Rrs_inf, ars1, ars2, ko, aw, aP_star, wavelength, bbw, bb_x)

def simulate_rrs_pixel(params, theta_s, theta_v, Rrs_inf, ars1, ars2, ko, aw, aP_star, wavelength, bbw, bb_x):
    import numpy as np
    zB, CX, CP, aY0, R_B = params
    bb = bbw + bb_x * CX
    a = aw + aP_star * CP + aY0 * np.exp(-0.011 * (wavelength - 440))
    wb = bb / (a + bb)
    Kd = ko * (a + bb) / np.cos(theta_s)
    Ku_W = (a + bb) * (1 + wb) ** k1u * (1 + k2u / np.cos(theta_s)) / np.cos(theta_v)
    Ku_B = (a + bb) * (1 + wb) ** k1b * (1 + k2b / np.cos(theta_s)) / np.cos(theta_v)
    exp_KuW = np.exp(-(Kd + Ku_W) * zB)
    exp_KuB = np.exp(-(Kd + Ku_B) * zB)
    Rrs_sim = Rrs_inf * (1 - ars1 * exp_KuW) + ars2 * (R_B / np.pi) * exp_KuB
    return Rrs_sim

def residual_zB(zB, band_values):
    return np.array([zB[0] - val for val in band_values])

def residual_CX(CX, band_values):
    return np.array([CX[0] - val for val in band_values])

def residual_CP(CP, band_values):
    return np.array([CP[0] - val for val in band_values])

def residual_aY(aY, band_values):
    return np.array([aY[0] - val for val in band_values])