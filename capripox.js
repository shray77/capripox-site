/*!
 * capripox.js — движок витрины «Оспа овец и коз: сценарии вакцинации».
 * DOM-free: чистая математика и данные, ноль обращений к документу —
 * гейтится known-answer'ами (tests/capripox.test.mjs, CI на каждый пуш).
 *
 * ЕДИНЫЙ ИСТОЧНИК всех чисел: EFSA AHAW Panel (2026). Sheep and goat pox:
 * vaccines and vaccination scenarios to control the epidemics in Greece and
 * Bulgaria in 2025. EFSA Journal 2026;24(6):e9928. DOI 10.2903/j.efsa.2026.9928
 * (открытый доступ; PMC12936713). Дополнительно: презентация WOAH RR-Europe
 * «SGE PPR and SGP, Mar 2026» (A. Broglia, EFSA BIOHAW) — 12 стр., параметры
 * модели слайдами. Каждая константа несёт поле src — секция/таблица статьи.
 */

/* ---------- 1. Ядро передачи (стохастическая пространственная модель BG) ---- */
/* Сила инфекции фермы j:  λ_j(t) = h · S_j(t) · Σ_{k≠j} I_k(t) · K(x_jk)
 * Ядро:                    K(x) = (1 + x/d0)^(−α)
 * Параметры — MLE по вспышкам Болгарии после 29.06.2025 в радиусе 50 км от
 * интродукции (статья, §2.2.2): h = 4.73·10⁻⁶, d0 = 1.25 км, α = 2.79.     */
export const KERNEL_BG = Object.freeze({
  h: 4.73e-6, d0: 1.25, alpha: 2.79,
  src: 'EFSA 2026.9928 §2.2.2 (метод максимального правдоподобия, Болгария)',
});

/** K(x) — относительный риск передачи на расстоянии x (км). K(0)=1. */
export function kernelK(x, p = KERNEL_BG) {
  if (!Number.isFinite(x) || x < 0 || !Number.isFinite(p.d0) || p.d0 <= 0 ||
      !Number.isFinite(p.alpha) || p.alpha <= 2) return NaN;
  return Math.pow(1 + x / p.d0, -p.alpha);
}

/** Накопленная доля силы инфекции от ферм в радиусе R (км) при равномерной
 *  плотности ферм: S(R) = ∫₀^R 2πx·K(x)dx / ∫₀^∞ 2πx·K(x)dx.
 *  Замкнутая форма (α>2): ∫x·K dx = d0²·[u^{2−α}/(2−α) − u^{1−α}/(1−α)], u=1+x/d0.
 *  Эталон сверен с численным интегрированием (scripts/capripox_reference.py,
 *  расхождение < 1e-13). Оговорка витрины: плотность реальная НЕ равномерна,
 *  MLE-оценка ядра исключала прыжки >50 км — профиль показываем как геометрию
 *  ядра, а не как готовый прогноз «сколько заражений случится в кольце».     */
export function ringShare(R, p = KERNEL_BG) {
  if (!Number.isFinite(R) || R < 0) return NaN;
  const { d0, alpha } = p;
  if (!(d0 > 0) || !(alpha > 2)) return NaN;
  const F = (u) => d0 * d0 * (Math.pow(u, 2 - alpha) / (2 - alpha) - Math.pow(u, 1 - alpha) / (1 - alpha));
  // Тотал = −F(1) = d0²·(1/(α−2) − 1/(α−1)). ГРАБЛЯ, пойманная гейтом
  // (S(1e9)→1): в первой версии между членами был ПЛЮС — тотал завышен
  // в 2,58 раза, и все доли колец занижены во столько же. Питон-эталон
  // подтверждал ошибку циркулярно (оба метода делились на один тотал);
  // чистая сверка — Симпсон + аналитический хвост — сходится с формой
  // до 1e-14 (scripts/capripox_reference.py, ред. 2).
  const total = d0 * d0 * (1 / (alpha - 2) - 1 / (alpha - 1));
  const part = F(1 + R / d0) - F(1);
  const s = part / total;
  return Math.min(1, Math.max(0, s)); // численная гигиена границ
}

/** Профиль ядра по набору радиусов (км): [{r, K, share}] — данные для таблицы. */
export function kernelProfile(radii, p = KERNEL_BG) {
  return radii.map((r) => ({ r, K: kernelK(r, p), share: ringShare(r, p) }));
}

/* ---------- 2. Финальный размер вспышки в стаде ----------------------------- */
/* Статья, §2.2.2: число заболевших в стаде — биномиальная выборка с вероятностью
 * f — финальным размером, корнем уравнения 1 − f − exp(−f·R0) = 0
 * (SIR, Keeling & Rohani 2008). Внутристадный R0 в статье НЕ задан — слайдер
 * пользователя; график честно помечен «R0 не задан в EFSA (2026)».            */
export function finalSize(R0) {
  if (!Number.isFinite(R0) || R0 <= 1) return NaN;
  const g = (f) => 1 - f - Math.exp(-f * R0);
  let lo = 0, hi = 1;
  for (let i = 0; i < 200; i++) {
    const mid = (lo + hi) / 2;
    if (g(mid) > 0) lo = mid; else hi = mid;
  }
  return (lo + hi) / 2;
}

/* ---------- 3. Пороговые покрытия и логистика кампании --------------------- */
/* EFSA AHAW (2014): контроль вспышек эффективен при покрытии ≥75% (и удержании).
 * Модель EFSA-2026 берёт VE = 95% (EFSA 2018), иммунитет 1 год, старт +90 дней
 * после интродукции, завершение вакцинации всех ферм за 60 дней (оценка BG-служб),
 * ежегодная замена поголовья 20% (оценка национальных служб).                 */

/** Минимальное покрытие для доли защищённых `target` при эффективности `ve`. */
export function minCoverage(target, ve) {
  if (!Number.isFinite(target) || !Number.isFinite(ve) || ve <= 0 || ve > 1 ||
      target <= 0 || target > ve) return null;
  return target / ve;
}

/** Дозы на голову: живые аттенуированные — 1 доза/животное (допущение модели
 *  EFSA: вакцинируются все животные фермы). Возвращает null на мусор.        */
export function campaign({ heads, coverage, ve = 0.95, farms = null }) {
  const h = Number(heads), c = Number(coverage), v = Number(ve);
  if (!Number.isFinite(h) || h < 0) return null;
  if (!Number.isFinite(c) || c <= 0 || c > 1) return null; // покрытие 0 — не кампания
  if (!Number.isFinite(v) || v < 0 || v > 1) return null;
  const out = {
    doses: Math.round(h * c),
    protectedShare: c * v,
    protectedHeads: Math.round(h * c * v),
    hit75: c * v >= 0.75,
    minCoverageFor75: minCoverage(0.75, v),
  };
  if (Number.isFinite(farms) && farms >= 0) out.farmsToVaccinate = Math.round(farms * c);
  return out;
}

/** Возврат восприимчивости после остановки кампании (модель EFSA: иммунитет
 *  вакцинированных 1 год; переболевшие (recovered) — пожизненно; замена
 *  поголовья 20%/год доезджает новыми восприимчивыми). Месяцы 0–24.
 *  immuneShare(months) = recovered + cov·ve·[months ≤ 12]                     */
export function susceptibleAfterStop(months, { coverage, ve = 0.95, recovered = 0 }) {
  if (!Number.isFinite(months) || months < 0) return null;
  const c = Number(coverage), r = Number(recovered);
  if (!Number.isFinite(c) || c < 0 || c > 1 || !Number.isFinite(r) || r < 0 || r > 1) return null;
  const immune = Math.min(1, r + c * ve * (months <= 12 ? 1 : 0));
  return 1 - immune;
}

/* ---------- 4. Данные: эпидемия, популяция, вакцины (EFSA 2026.9928) ------- */

export const EPIDEMIC = Object.freeze({
  window: '1 января 2024 — 29 января 2026 (ADIS)',
  countries: [
    { c: 'Греция',      outbreaks: 2081, note: 'мультикластерная эпидемия, средний межкластерный разброс ~200 км', src: '§3.2.2' },
    { c: 'Болгария',    outbreaks: 203,  note: 'единичная интродукция, распространение ~180 км от границы', src: '§3.2.2' },
    { c: 'Румыния',     outbreaks: 27,   note: 'с июня 2025', src: '§3.2.2' },
    { c: 'Сербия',      outbreaks: 3,    note: '2 в сентябре + 1 в октябре 2025', src: '§3.2.2' },
    { c: 'Сев. Македония', outbreaks: 1, note: '23.01.2026 (подтверждена 27.01)', src: '§3.2.2' },
  ],
  year2024: { bg: 12, gr: 328, total: 340, src: '§3.2.2 (ADIS)' },
  peak2025: 354,           /* пик месячного числа вспышек в Греции, 2025 */
  bg2025byRegion: [
    { r: 'Пловдив', n: 144 }, { r: 'Хасково', n: 16 }, { r: 'Пазарджик', n: 15 },
    { r: 'Стара-Загора', n: 10 }, { r: 'Сливен', n: 3 }, { r: 'Бургас', n: 1 },
    { r: 'Кырджали', n: 1 }, { r: 'Кюстендил', n: 1 },
  ], /* на 05.12.2025, 8 регионов, 191 вспышка — §3.2.2 */
  clusters: {
    method: 'ST-DBSCAN: ≤28 дней, ≤20 км; минимум 25 точек (GR) / 5 (BG, RO)',
    greece: 10, bulgaria: 6, romania: 1, src: '§3.2.3, Table 3',
  },
  origin: 'восток на запад от границы Турции; штамм почти идентичен во всех странах (NRL GR/BG/RO + EURL)',
});

export const POPULATION = Object.freeze({
  greece:  { farms: 83500, heads: 17000000, sheep: 13000000, goats: 4000000,
             densityHeadsKm2: 129, densityFarmsKm2: 0.64, src: '§3.2.1' },
  bulgaria:{ farms: 24300, heads: 1400000, sheep: 1230000, goats: 200000,
             densityHeadsKm2: 13, densityFarmsKm2: 0.22, src: '§3.2.1' },
  zones: [
    { id: 'bg10', label: 'Болгария, зоны наблюдения 10 км',  farms: 2500,  heads: 192000, src: 'Table 4' },
    { id: 'bg20', label: 'Болгария, зоны наблюдения 20 км',  farms: 5000,  heads: 370000, src: 'Table 4' },
    { id: 'bg4r', label: 'Болгария, 4 региона (Пловдив, Хасково, Стара-Загора, Пазарджик)',
      farms: 4078, heads: 258953, src: 'Table 5' },
    { id: 'gr20', label: 'Греция, 20-км зоны (все затронутые РО)', farms: 30392, heads: 5613183, src: 'Table 6 (36% ферм, 33% голов)' },
    { id: 'grA',  label: 'Греция, затронутые регионы (к 17.12.2025)', farms: 54000, heads: 9900000, src: '§3.2.4.1' },
    { id: 'grW1', label: 'Греция, западный приоритет — регионы целиком', farms: null, heads: 4400000, src: '§3.2.4.2' },
    { id: 'grW2', label: 'Греция, западный приоритет — 20-км зоны', farms: null, heads: 2500000, src: '§3.2.4.2' },
  ],
});

export const VACCINES = Object.freeze({
  litReview: {
    scope: 'Elicit-обзор: 497 → 139 бумаг (эффективность), 174 → 113 (безопасность); извлечение по 112+27 и 54 полным текстам',
    liveAttenuated: 'эффективность 80–100% в большинстве исследований; иммунитет ≥12 мес у овец, до 52 мес у коз',
    inactivated: 'до 9 мес (Boumart 2016); полная защита в части исследований на высоких дозах',
    maternal: '~60 дней у ягнят (Ebrahimi-Jam 2022)',
    noDIVA: 'различающих вакцин (DIVA) на рынке нет; в ЕС авторизованных вакцин нет',
    src: '§3.1.1, §3.1.1.2, §3.1.1.3',
  },
  eurlSafety: { /* Table 1: испытания EURL, 2 коммерческие вакцины на штамм */
    head: ['Штамм', 'Фебр. >41°C', 'Тяж. местн. реакции', 'Назальн. шеддинг', 'Вирусемия'],
    rows: [
      ['RM-65',    '26/36 (72%)', '18/36 (50%)', '7/36 (19%), ~1 д', '8/36 (22%)'],
      ['Romania',  '23/28 (82%)', '17/28 (61%)', '4/28 (14%), ~1,5 д', '10/28 (36%)'],
      ['Bakirköy', '1/28 (4%)',   '0/28 (0%)',   '17/28 (61%), ~1,1 д', '11/28 (39%)'],
    ],
    note: 'доза >10²·⁵ TCID₅₀ (минимум WOAH); шеддинг/вирусемия короткие и близки к пределу детекции',
    src: '§3.1.2, Table 1',
  },
  eurlEfficacy: { /* Table 2: челлендж вирулентным марокканским SPPV, 10^3.3–10^4.95 TCID₅₀, интра-дермально, 21 день */
    head: ['Группа', 'Смертность', 'Морбидность', 'Шеддинг челленджа'],
    rows: [
      ['Контроль (n=20)', '20 (100%)', '20 (100%)', 'высокий'],
      ['RM-65 (n=18)',    '0 (0%)',    '0 (0%)',    '~1000-кратное снижение'],
      ['Romania (n=18)',  '0 (0%)',    '0 (0%)',    '~1000-кратное снижение'],
      ['Bakirköy (n=28)', '0 (0%)',    '6 (21%) лёгкая', '~1000-кратное снижение'],
    ],
    note: 'жёсткий челлендж (высоковирулентный штамм, высокая доза, интра-дермально) — полевые условия мягче',
    src: '§3.1.2, Table 2',
  },
  keyFacts: [
    'все три штамма — 100% защита от смертности; RM-65 и Romania — и от morbidity',
    'Bakirköy — самый безопасный (фебр. 4% против 72–82%), ценой лёгкой morbidity 21% при жёстком челлендже',
    'риск передачи вакцинного вируса минимален; держать животных на ферме 3 недели после вакцинации',
    'живые вакцины RM-65/Romania могут давать лихорадку и местные реакции — предупреждать фермеров заранее',
  ],
});

export const MODEL = Object.freeze({
  type: 'стохастическая пространственно-ядровая, S-I-R-V на уровне ферм',
  durationScenarios: [
    { d: 15, meaning: 'быстрая детекция и зачистка: инкубация 10 д + 5 д от репорта до зачистки' },
    { d: 30, meaning: 'нормальная детекция: 25 д (2–3 инкубационных периода) + 5 д' },
    { d: 60, meaning: 'без зачистки: средняя длительность внутристадной вспышки' },
  ],
  vaccinationAssumptions: {
    ve: 0.95,                /* EFSA 2018 */
    immunityYears: 1,        /* остановка кампании → вакцинированные снова восприимчивы */
    startDay: 90,            /* старт кампании после интродукции вируса */
    completeDays: 60,        /* постоянный темп, все фермы привиты за 60 дней (оценка BG) */
    replacement: 0.20,       /* ежегодная замена поголовья (оценка нацслужб) */
    scope: '4 региона (Пловдив, Хасково, Стара-Загора, Пазарджик) или вся страна',
  },
  greeceCaveat: 'для Греции ядро не удалось оценить (мультикластерность, горы/острова) — выводы по Греции описательные, выводы Болгарии применимы как минимум',
});

/* ---------- 5. Матрица стратегий (выводы §3.2.4.2 + Conclusions) ------------ */
/* Ключ: det ∈ {15,30,60} дней до зачистки; vac ∈ none|regions4|country;
 * ctrl ∈ ''|'20km' (жёсткий 20-км контроль перемещений). Проверяются только
 * комбинации, которые реально описаны в статье; прочие честно помечены
 * modeled:false — никаких экстраполяций от себя.                              */

const _S = {
  det30_none: 'Быстрая детекция и зачистка (~30 дней после интродукции: 25 д + 5 д от репорта до зачистки), проводимая систематически, — достаточна, чтобы взять эпидемию под контроль за 1–2 года.',
  det60_none: 'Без контроля модель предсказала распространение SGP по всей Болгарии и переход в эндемию.',
  det30_country: 'Общенациональная вакцинация значимо снижала и число инфицированных ферм, и географию распространения (строже, чем детекция+зачистка без вакцинации). Двухлетняя кампания на фоне нормальной детекции/зачистки достаточна для элиминации SGP — на 1 год раньше, чем без вакцинации.',
  det30_r4_20km: 'Вакцинация только 4 изначально затронутых регионов снимает риск выхода вируса из вакцинированной зоны ПРИ жёстком 20-км контроле перемещений внутри зон наблюдения (в модели — нулевая передача далее 20 км от инфицированной фермы).',
  det30_r4_no: 'Вакцинация только 4 регионов без жёсткого 20-км контроля — есть риск выхода вируса за пределы вакцинированной области и более широкого распространения.',
  det60_country: 'Если других мер нет, вакцинация становится определяющей для контроля: нужна кампания не менее 3 лет; элиминация достигается за 2–3 года. Вывод EFSA: оптимальна комбинация с ранней детекцией, быстрой зачисткой и контролем перемещений.',
  notModeled_det15: 'Сценарий 15 дней закладывался в модель как «быстрая детекция» (1 инкубационный период + 5 д), но отдельного вывода EFSA для него нет — выводы сформулированы для ~30 дней («в пределах ~30 дней систематически»).',
  notModeled_combo: 'Такое сочетание в EFSA (2026) не моделировалось. Статья описывает региональную вакцинацию на фоне стандартных мер детекции/зачистки — перенос на другие комбинации был бы экстраполяцией.',
};

export function scenarioVerdict(det, vac, ctrl) {
  const d = Number(det), c20 = ctrl === '20km';
  if (![15, 30, 60].includes(d)) return null;
  if (!['none', 'regions4', 'country'].includes(vac)) return null;
  if (vac === 'none') {
    if (d === 15) return { modeled: false, verdict: _S.notModeled_det15, src: '§2.2.2, Conclusions' };
    if (d === 30) return { modeled: true, verdict: _S.det30_none, src: '§3.2.4.2, Conclusions' };
    return { modeled: true, verdict: _S.det60_none, src: '§3.2.4.2' };
  }
  if (d === 60 && vac === 'country') {
    /* «Если других мер нет» — сценарий статьи: вакцинация без детекции/зачистки */
    return { modeled: true, verdict: _S.det60_country, src: '§3.2.4.2, Conclusions' };
  }
  if (d !== 30) {
    return { modeled: false, verdict: _S.notModeled_combo, src: '§3.2.4.2 (охват модельных сценариев)' };
  }
  if (vac === 'country') return { modeled: true, verdict: _S.det30_country, src: '§3.2.4.2, Conclusions' };
  return c20
    ? { modeled: true, verdict: _S.det30_r4_20km, src: '§3.2.4.2, Conclusions' }
    : { modeled: true, verdict: _S.det30_r4_no, src: '§3.2.4.2, Conclusions' };
}

/** Перечень всех легитимных ключей матрицы (для гейта полноты в тестах). */
export function scenarioKeys() {
  const out = [];
  for (const d of [15, 30, 60])
    for (const v of ['none', 'regions4', 'country'])
      for (const c of ['', '20km'])
        out.push([d, v, c]);
  return out;
}

/* ---------- 6. Разбор пользовательского ввода (уроки Task 178) ------------- */

/** Число голов «своего» региона: запятая/точка, пробелы-разделители, мусор
 *  перечисляется явно. Отрицательные и не-числа — ошибка.                     */
export function parseHeads(input) {
  if (typeof input === 'number') {
    if (!Number.isFinite(input) || input < 0) return { ok: false, bad: [String(input)] };
    return { ok: true, value: input };
  }
  const s = String(input ?? '').trim();
  if (s === '') return { ok: false, bad: ['пусто'] };
  const norm = s.replace(/\s+/g, '').replace(',', '.');
  const v = Number(norm);
  if (!Number.isFinite(v) || v < 0) return { ok: false, bad: [s.slice(0, 30)] };
  if (v > 1e9) return { ok: false, bad: [s.slice(0, 30) + ' — больше 10⁹ голов?'] };
  return { ok: true, value: v };
}

export const SOURCES = Object.freeze([
  { t: 'EFSA AHAW Panel (2026). Sheep and goat pox: vaccines and vaccination scenarios to control the epidemics in Greece and Bulgaria in 2025. EFSA Journal 24(6):e9928', u: 'https://doi.org/10.2903/j.efsa.2026.9928', tag: 'первоисточник, открытый доступ' },
  { t: 'Тот же текст в Europe PMC (PMC12936713, PMID 41768356)', u: 'https://europepmc.org/article/MED/41768356', tag: 'зеркало' },
  { t: 'Broglia A. (EFSA BIOHAW). SGE PPR and SGP, Mar 2026: Sheep and goat pox — vaccines and vaccination strategies (WOAH RR-Europe, 12 стр.)', u: 'https://rr-europe.woah.org/app/uploads/2026/04/6_SPGP_vaccines_EFSA.pdf', tag: 'презентация параметров модели' },
  { t: 'EFSA. Sheep pox and goat pox — Disease profile (карточка болезни: летальность по группам)', u: 'https://animal-diseases.efsa.europa.eu/SPPV', tag: 'профайл' },
  { t: 'Приложение EFSA: лит-обзор эффективности вакцин (Annex 1, Zenodo)', u: 'https://www.zenodo.org/uploads/18439931', tag: 'приложение' },
  { t: 'Приложение EFSA: лит-обзор безопасности (Annex 2, Zenodo)', u: 'https://www.zenodo.org/uploads/18374862', tag: 'приложение' },
  { t: 'Динамическая карта распространения SGP 2024–2026 (EFSA, Zenodo)', u: 'https://www.zenodo.org/uploads/18439735', tag: 'анимация' },
  { t: 'EFSA AHAW Panel (2014). Scientific opinion on sheep and goat pox (порог покрытия 75%)', u: 'https://doi.org/10.2903/j.efsa.2014.3885', tag: 'порог 75%' },
  { t: 'Regulation (EU) 2020/687 — правила профилактики и контроля занесённых болезней (зоны наблюдения 10→20 км)', u: 'https://eur-lex.europa.eu/legal-content/EN/TXT/?uri=CELEX:32020R0687', tag: 'право ЕС' },
]);
