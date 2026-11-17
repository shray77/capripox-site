/*!
 * capripox.test.mjs — гейт витрины «Оспа овец и коз»: известные ответы.
 *
 * Классы проверок (гейт-культура экосистемы, паттерн fecrt.test.mjs):
 *  1. ядро K(x) — значения, замороженные scripts/capripox_reference.py
 *     (Python-эталон, сам сверенный двумя методами: закрытая форма и
 *     численное интегрирование, расхождение <1e-13);
 *  2. кольца S(R) — арифметика против того же эталона + монотонность;
 *  3. финальный размер — корни 1−f−exp(−f·R0)=0 (бисекция vs сканирование);
 *  4. логистика кампании — точечная арифметика (дозы, порог 75%, VE 95%);
 *  5. матрица стратегий — полнота (12 комбинаций) + точные подстроки
 *     выводов EFSA + честность немоделированных клеток;
 *  6. данные — числа статьи (203/2081/27/3/1, популяции, EURL) против
 *     референса tests/reference.json — дрейф правок = красный CI;
 *  7. разбор ввода — мусор, запятые, края.
 *
 * Запуск: node --test tests/capripox.test.mjs   (CI: .github/workflows/ci.yml)
 */
import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { dirname, join } from 'node:path';

import * as C from '../capripox.js';

const HERE = dirname(fileURLToPath(import.meta.url));
const near = (a, b, eps = 1e-9) => Math.abs(a - b) < eps;

/* ---------- 1. ядро K(x): замороженный Python-эталон ----------------------- */

test('K(x): эталонные точки (d0=1.25, α=2.79, MLE Болгарии)', () => {
  assert.ok(near(C.kernelK(0), 1.0));
  assert.ok(near(C.kernelK(0.5), 3.911134763834481e-1));
  assert.ok(near(C.kernelK(1), 1.939946310936210e-1));
  assert.ok(near(C.kernelK(1.25), 1.445860229881609e-1));   // 2^(−2.79)
  assert.ok(near(C.kernelK(2.5), 4.664785724232973e-2));
  assert.ok(near(C.kernelK(5), 1.121692169022263e-2));
  assert.ok(near(C.kernelK(10), 2.176022585300775e-3));
  assert.ok(near(C.kernelK(20), 3.690179648358441e-4));
  assert.ok(near(C.kernelK(50), 3.164703590327204e-5));
  assert.ok(near(C.kernelK(100), 4.735074291739067e-6));
});

test('K(x): монотонное убывание и мусор на входе', () => {
  let prev = C.kernelK(0);
  for (let x = 1; x <= 60; x += 1) {
    const k = C.kernelK(x);
    assert.ok(k < prev, `K не убывает в x=${x}`);
    prev = k;
  }
  assert.ok(Number.isNaN(C.kernelK(-1)));
  assert.ok(Number.isNaN(C.kernelK(NaN)));
  assert.ok(Number.isNaN(C.kernelK(5, { d0: 1.25, alpha: 1.5 }))); // α≤2 — не интегрируемо
  assert.ok(Number.isNaN(C.kernelK(5, { d0: 0, alpha: 2.79 })));
});

test('K(x): чувствительность к параметрам (слайдеры витрины не ломают математику)', () => {
  assert.ok(near(C.kernelK(20, { d0: 1.0, alpha: 2.79 }), 2.046483e-04, 1e-8));
  assert.ok(near(C.kernelK(20, { d0: 1.25, alpha: 3.1 }), 1.533234e-04, 1e-8));
  assert.ok(C.kernelK(20, { d0: 1.5, alpha: 2.79 }) > C.kernelK(20)); // шире ядро → дальше риск
});

/* ---------- 2. кольца S(R): доля силы инфекции в радиусе ------------------- */

test('S(R): эталон (равномерная плотность ферм, оговорка в витрине)', () => {
  assert.ok(near(C.ringShare(0), 0));
  assert.ok(near(C.ringShare(0.5), 6.0388984336404e-2, 1e-12));
  assert.ok(near(C.ringShare(1), 1.50769102924565e-1, 1e-12));
  assert.ok(near(C.ringShare(2), 3.01389956602315e-1, 1e-12));
  assert.ok(near(C.ringShare(3), 4.07619553249438e-1, 1e-12));
  assert.ok(near(C.ringShare(5), 5.42349595038917e-1, 1e-12));
  assert.ok(near(C.ringShare(10), 6.99970005938729e-1, 1e-12));
  assert.ok(near(C.ringShare(20), 8.14059227878515e-1, 1e-12));
  assert.ok(near(C.ringShare(50), 9.05799432930321e-1, 1e-12));
  assert.ok(near(C.ringShare(100), 9.44693385257629e-1, 1e-12));
});

test('S(R): монотонность, границы, мусор', () => {
  let prev = -1;
  for (let R = 0; R <= 80; R += 2) {
    const s = C.ringShare(R);
    assert.ok(s >= prev - 1e-12 && s <= 1, `S(R) не монотонна/выходила за 1 при R=${R} (допуск FP 1e-12)`);
    prev = s;
  }
  assert.ok(near(C.ringShare(1e9), 1, 1e-6)); // хвост интегрируем при α>2
  assert.ok(Number.isNaN(C.ringShare(-5)));
  assert.ok(Number.isNaN(C.ringShare(NaN)));
});

test('kernelProfile: радиусы → K и share согласованы', () => {
  const p = C.kernelProfile([3, 10, 20]);
  assert.deepEqual(p.map((x) => x.r), [3, 10, 20]);
  assert.ok(near(p[2].K, C.kernelK(20)));
  assert.ok(near(p[2].share, C.ringShare(20)));
});

/* ---------- 3. финальный размер стада -------------------------------------- */

test('finalSize: корни 1−f−exp(−f·R0)=0 (эталон бисекции+сканирования)', () => {
  assert.ok(near(C.finalSize(1.05), 0.093701837072901, 1e-8));
  assert.ok(near(C.finalSize(1.2), 0.313698331041217, 1e-8));
  assert.ok(near(C.finalSize(1.5), 0.582811643865811, 1e-8));
  assert.ok(near(C.finalSize(2), 0.796812130020020, 1e-8));
  assert.ok(near(C.finalSize(3), 0.940479790707360, 1e-8));
  assert.ok(near(C.finalSize(4), 0.980172598718221, 1e-8));
  assert.ok(near(C.finalSize(5), 0.993022846348855, 1e-8));
  assert.ok(near(C.finalSize(8), 0.999663633449189, 1e-8));
  assert.ok(near(C.finalSize(10), 0.999954579444654, 1e-8));
  assert.ok(near(C.finalSize(20), 0.999999997938846, 1e-8));
});

test('finalSize: подстановка в уравнение даёт ноль; R0≤1 — NaN (нет эпидемии)', () => {
  for (const R0 of [1.5, 2, 3.7, 5, 12]) {
    const f = C.finalSize(R0);
    assert.ok(Math.abs(1 - f - Math.exp(-f * R0)) < 1e-10, `невязка при R0=${R0}`);
  }
  assert.ok(Number.isNaN(C.finalSize(1)));
  assert.ok(Number.isNaN(C.finalSize(0.5)));
  assert.ok(Number.isNaN(C.finalSize(NaN)));
});

/* ---------- 4. пороги и логистика кампании -------------------------------- */

test('minCoverage: порог 75% защищённых при VE 95% — 78,9% покрытия', () => {
  assert.ok(near(C.minCoverage(0.75, 0.95), 0.75 / 0.95));
  assert.ok(near(C.minCoverage(0.75, 0.95), 0.7894736842105263));
  assert.equal(C.minCoverage(0.75, 0.95) > 0.75, true); // потери на VE
  assert.ok(near(C.minCoverage(0.8, 0.95), 0.8421052631578948)); // 80% защищённых — 84,2% покрытия
  assert.equal(C.minCoverage(0.96, 0.95), null);         // цель выше VE — недостижимо
  assert.equal(C.minCoverage(0, 0.95), null);
});

test('campaign: дозы и защищённая доля (арифметика руками)', () => {
  const r = C.campaign({ heads: 258953, coverage: 0.8, ve: 0.95, farms: 4078 });
  assert.equal(r.doses, 207162);            // 258953×0.8 = 207162.4 → округление
  assert.ok(near(r.protectedShare, 0.76));
  assert.equal(r.protectedHeads, 196804);   // 258953×0.8×0.95 = 196804.28
  assert.equal(r.hit75, true);
  assert.equal(r.farmsToVaccinate, 3262);   // 4078×0.8
  const r2 = C.campaign({ heads: 5613183, coverage: 0.75, ve: 0.95 });
  assert.equal(r2.doses, 4209887);          // 5613183×0.75 = 4209887.25
  assert.equal(r2.hit75, false);            // 0.75×0.95 = 71.25% < 75%
  assert.equal(r2.minCoverageFor75, C.minCoverage(0.75, 0.95));
});

test('campaign: мусор и края', () => {
  assert.equal(C.campaign({ heads: 'abc', coverage: 0.8 }), null);
  assert.equal(C.campaign({ heads: -100, coverage: 0.8 }), null);
  assert.equal(C.campaign({ heads: 1000, coverage: 1.5 }), null);
  assert.equal(C.campaign({ heads: 1000, coverage: 0 }), null);
  assert.equal(C.campaign({ heads: 0, coverage: 0.8 }).doses, 0);
  const full = C.campaign({ heads: 1000, coverage: 1, ve: 0.95 });
  assert.equal(full.doses, 1000);
  assert.ok(near(full.protectedShare, 0.95));
});

test('susceptibleAfterStop: иммунитет 1 год, потом восприимчивость возвращается', () => {
  const s0 = C.susceptibleAfterStop(0, { coverage: 0.8, ve: 0.95 });
  assert.ok(near(s0, 1 - 0.76));
  const s12 = C.susceptibleAfterStop(12, { coverage: 0.8, ve: 0.95 });
  assert.ok(near(s12, 1 - 0.76));           // ровно 12 мес — ещё защищены
  const s13 = C.susceptibleAfterStop(13, { coverage: 0.8, ve: 0.95 });
  assert.ok(near(s13, 1));                  // год прошёл — иммунитета нет
  const sR = C.susceptibleAfterStop(24, { coverage: 0.8, ve: 0.95, recovered: 0.1 });
  assert.ok(near(sR, 0.9));                 // переболевшие остаются иммунными
  assert.equal(C.susceptibleAfterStop(-1, { coverage: 0.8 }), null);
  assert.equal(C.susceptibleAfterStop(6, { coverage: 2 }), null);
});

/* ---------- 5. матрица стратегий ------------------------------------------ */

test('матрица: полнота — у всех 12 комбинаций есть вердикт и источник', () => {
  for (const [d, v, c] of C.scenarioKeys()) {
    const r = C.scenarioVerdict(d, v, c);
    assert.ok(r, `нет вердикта для ${d}/${v}/${c}`);
    assert.ok(typeof r.verdict === 'string' && r.verdict.length > 20);
    assert.ok(typeof r.src === 'string' && r.src.includes('EFSA') === false ? r.src.length > 3 : true);
    assert.equal(typeof r.modeled, 'boolean');
  }
});

test('матрица: ключевые клетки — точные формулировки выводов EFSA', () => {
  const a = C.scenarioVerdict(30, 'none', '');
  assert.equal(a.modeled, true);
  assert.ok(a.verdict.includes('1–2 года'));
  const b = C.scenarioVerdict(60, 'none', '');
  assert.equal(b.modeled, true);
  assert.ok(b.verdict.includes('эндемию'));
  const c = C.scenarioVerdict(30, 'country', '');
  assert.equal(c.modeled, true);
  assert.ok(c.verdict.includes('1 год раньше'));
  const d = C.scenarioVerdict(30, 'regions4', '20km');
  assert.equal(d.modeled, true);
  assert.ok(d.verdict.includes('20-км'));
  const e = C.scenarioVerdict(30, 'regions4', '');
  assert.equal(e.modeled, true);
  assert.ok(e.verdict.includes('риск выхода'));
  const f = C.scenarioVerdict(60, 'country', '');
  assert.equal(f.modeled, true);
  assert.ok(f.verdict.includes('3 лет') && f.verdict.includes('2–3 года'));
});

test('матрица: немоделированные клетки помечены честно, без отсебятины', () => {
  for (const [d, v, c] of [[15, 'none', ''], [15, 'none', '20km'], [15, 'country', ''],
                           [60, 'regions4', ''], [60, 'regions4', '20km']]) {
    const r = C.scenarioVerdict(d, v, c);
    assert.equal(r.modeled, false, `${d}/${v}/${c} должен быть немоделированным`);
    assert.ok(r.verdict.includes('не моделировал') || r.verdict.includes('отдельного вывода'),
      `нет честной пометки для ${d}/${v}/${c}`);
  }
  assert.equal(C.scenarioVerdict(45, 'none', ''), null);   // несуществующий сценарий
  assert.equal(C.scenarioVerdict(30, 'moon', ''), null);
});

/* ---------- 6. данные статьи против референса ------------------------------ */

const REF = JSON.parse(readFileSync(join(HERE, 'reference.json'), 'utf8'));

test('данные: вспышки по странам — точные числа статьи', () => {
  const byC = Object.fromEntries(C.EPIDEMIC.countries.map((x) => [x.c, x.outbreaks]));
  assert.deepEqual(
    { 'Греция': byC['Греция'], 'Болгария': byC['Болгария'], 'Румыния': byC['Румыния'], 'Сербия': byC['Сербия'], 'Сев. Македония': byC['Сев. Македония'] },
    REF.outbreaksByCountry,
  );
  assert.equal(C.EPIDEMIC.year2024.total, 340);
  assert.equal(C.EPIDEMIC.peak2025, 354);
  const bgSum = C.EPIDEMIC.bg2025byRegion.reduce((s, x) => s + x.n, 0);
  assert.equal(bgSum, 191); // §3.2.2: 191 вспышка в 8 регионах на 05.12.2025
});

test('данные: популяции и зоны (Table 4–6)', () => {
  assert.equal(C.POPULATION.greece.farms, REF.greeceFarms);
  assert.equal(C.POPULATION.greece.heads, REF.greeceHeads);
  assert.equal(C.POPULATION.bulgaria.farms, REF.bulgariaFarms);
  assert.equal(C.POPULATION.bulgaria.heads, REF.bulgariaHeads);
  const z = Object.fromEntries(C.POPULATION.zones.map((x) => [x.id, x.heads]));
  assert.deepEqual(
    { bg10: z.bg10, bg20: z.bg20, bg4r: z.bg4r, gr20: z.gr20, grA: z.grA, grW1: z.grW1, grW2: z.grW2 },
    REF.zonesHeads,
  );
});

test('данные: параметры модели (MLE Болгарии) — сверка с формулой статьи', () => {
  assert.equal(C.KERNEL_BG.h, REF.h);
  assert.ok(near(C.KERNEL_BG.d0, REF.d0));
  assert.ok(near(C.KERNEL_BG.alpha, REF.alpha));
  assert.equal(C.MODEL.vaccinationAssumptions.ve, 0.95);
  assert.equal(C.MODEL.vaccinationAssumptions.immunityYears, 1);
  assert.equal(C.MODEL.vaccinationAssumptions.startDay, 90);
  assert.equal(C.MODEL.vaccinationAssumptions.completeDays, 60);
  assert.equal(C.MODEL.vaccinationAssumptions.replacement, 0.20);
});

test('данные: EURL-таблицы — контроль 100% смертности, Bakirköy 21% лёгкой морбидности', () => {
  const ctl = C.VACCINES.eurlEfficacy.rows.find((r) => r[0].startsWith('Контроль'));
  assert.equal(ctl[1], '20 (100%)');
  assert.equal(ctl[2], '20 (100%)');
  const bak = C.VACCINES.eurlEfficacy.rows.find((r) => r[0].startsWith('Bakirköy'));
  assert.ok(bak[2].includes('21%'));
  const safety = C.VACCINES.eurlSafety.rows;
  assert.equal(safety.length, 3);
  for (const row of safety) assert.equal(row.length, 5);
});

test('источники: все с URL, DOI первоисточника на месте', () => {
  assert.ok(C.SOURCES.length >= 9);
  for (const s of C.SOURCES) assert.ok(s.u.startsWith('https://'));
  assert.ok(C.SOURCES.some((s) => s.u.includes('10.2903/j.efsa.2026.9928')));
});

/* ---------- 7. разбор ввода: мусор, запятые, края -------------------------- */

test('parseHeads: число голов — RU-запятая, пробелы, тысячи', () => {
  assert.deepEqual(C.parseHeads('258953'), { ok: true, value: 258953 });
  assert.deepEqual(C.parseHeads('258 953'), { ok: true, value: 258953 });
  assert.deepEqual(C.parseHeads('25000,5'), { ok: true, value: 25000.5 });
  assert.deepEqual(C.parseHeads(1200), { ok: true, value: 1200 });
  assert.deepEqual(C.parseHeads('0'), { ok: true, value: 0 });
});

test('parseHeads: мусор и края падают явно, а не тихо', () => {
  assert.equal(C.parseHeads('abc').ok, false);
  assert.equal(C.parseHeads('-5').ok, false);
  assert.equal(C.parseHeads('').ok, false);
  assert.equal(C.parseHeads('1e12').ok, false);      // >10⁹ — явная ошибка
  assert.equal(C.parseHeads(NaN).ok, false);
  assert.equal(C.parseHeads(undefined).ok, false);
});
