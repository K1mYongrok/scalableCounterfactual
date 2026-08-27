# scalableCounterfactual 시작 안내

이 폴더는 특정 임금자료에 종속되지 않은 범용 R 패키지 소스입니다. 도시·농촌 임금격차 분석에 사용한 변수 생성, 표본 제한, 데이터 파일은 포함하지 않습니다.

## 가장 먼저 볼 파일

1. `README.md`: 패키지 개요, 지원 모형과 전체 사용법
2. `inst/examples/quick_start.R`: 데이터 생성부터 분해·요약·그림까지 실행되는 예제
3. `inst/doc/ARCHITECTURE.md`: 추정, marginalization, bootstrap 모듈 구조
4. `inst/doc/RELEASE_VALIDATION_1.0.0.md`: 현재 버전의 검증 범위와 남은 과제

## 설치

GitHub의 공개 소스를 설치하려면 R 콘솔에서 다음을 실행합니다.

```r
install.packages("remotes")
remotes::install_github("K1mYongrok/scalableCounterfactual")
library(scalableCounterfactual)
```

의존 패키지가 없으면 먼저 설치합니다.

```r
install.packages(c("data.table", "Hmisc", "quantreg", "survival"))
```

개발 중인 소스를 직접 설치하려면 이 폴더의 상위 경로에서 다음을 실행합니다.

```powershell
R CMD INSTALL scalableCounterfactual
```

CUDA는 선택 사항입니다. NVIDIA GPU가 없거나 CUDA 환경을 설치하지 않아도 CPU 기능은 모두 사용할 수 있습니다.

## 2분 실행 예제

```r
library(scalableCounterfactual)

set.seed(42)
n <- 600
d <- data.frame(
  outcome = rnorm(n),
  x1 = rnorm(n),
  x2 = rbinom(n, 1, 0.4),
  group = rep(0:1, each = n / 2),
  sampling_weight = runif(n, 0.5, 2)
)

fit <- counterfactual_decompose(
  outcome ~ x1 + x2,
  data = d,
  group = "group",
  weights = "sampling_weight",
  model = "qr",
  solver = "pfnb",
  control = cf_control(
    nreg = 9,
    reported_quantiles = c(0.25, 0.5, 0.75),
    bootstrap_progress = FALSE
  ),
  bootstrap_reps = 3,
  seed = 42
)

summary(fit)
plot(fit)
```

이 예제는 패키지 작동만 확인하는 범용 예제입니다. 실제 논문의 표본 제한,
변수 생성, 도시·농촌 정의와 회귀식은 패키지 기본값이 아니며 별도의 분석
스크립트에서 명시해야 합니다.

전체 예제는 다음 명령으로 실행할 수 있습니다.

```r
source(system.file("examples", "quick_start.R", package = "scalableCounterfactual"))
```

## 핵심 인수

| 인수 | 의미 | 대표 선택지 |
|---|---|---|
| `model` | 조건부분포 모형 | `qr`, `cqr`, `loc`, `locsca`, `logit`, `probit`, `cloglog`, `lpm`, `cox` |
| `solver` | QR 계산 알고리즘 | `br`, `fn`, `pfn`, `qfnb`, `pfnb`, `proqreg`, `profn`, `onestep`, `auto` |
| `bootstrap_reps` | bootstrap 반복 수 | 시험 10–100, 최종 분석은 연구 설계에 맞게 확대 |
| `point_workers` | 두 집단 point fit 병렬 수 | 일반적으로 1–2 |
| `bootstrap_workers` | bootstrap draw 병렬 수 | RAM에 맞게 설정 |
| `gpu_backend` | 예측·marginalization 장치 | `cpu`, `cuda`, `auto` |

`pfnb`는 표준 QR 목적함수를 푸는 exact solver 계열이고, `onestep`은 Stata `qrprocess` 1.1.3의 one-step 알고리즘을 추적 가능한 형태로 옮긴 근사 process solver입니다. CUDA는 현재 PFNB fitting 자체가 아니라 대규모 예측과 counterfactual marginalization을 가속합니다.

## 폴더 구조

```text
scalableCounterfactual/
├─ R/                 패키지의 추정·분해·inference 소스
├─ man/               함수 도움말
├─ tests/             API, 수치 동등성, GPU, Stata parity 테스트
├─ inst/
│  ├─ doc/            현재 설계·API·검증 기술 문서
│  ├─ examples/       설치 후 실행 가능한 예제
│  ├─ python/         선택적 CUDA backend와 환경 명세
│  ├─ provenance/     외부 알고리즘 출처와 버전 고정 정보
│  └─ scripts/        범용 명령행 실행 도구
├─ tools/             release 검증용 개발 스크립트
├─ DESCRIPTION        패키지 metadata와 의존성
├─ NAMESPACE          공개 API
├─ NEWS.md            버전별 변경사항
└─ README.md          전체 사용자 문서
```

## 결과 해석

기본 분해 방향은 패키지에 전달한 reference group과 comparison group에 의해 결정됩니다. 결과 객체와 `run_metadata.csv`에서 실제 집단 방향을 확인해야 합니다.

일반적인 출력은 다음을 포함합니다.

- `decomposition.csv`: total, structure, composition 효과와 신뢰구간
- `distribution_diagnostics.csv`: 관측·반사실 분포 진단
- `marginalization_diagnostics.csv`: CPU/CUDA 및 메모리 처리 경로
- `bootstrap_resources.csv`: draw별 시간, 실패 및 메모리 기록
- `functional_effect_tests.csv`: 구간별 KS/CMS 검정
- `run_metadata.csv`: 모형, solver, 장치, 표본 및 재현성 정보

## 현재 상태

현재 공유 버전은 `1.0.0`입니다. CPU와 실제 RTX 4060 Ti CUDA 검증, 분석적 simulation, solver benchmark 및 bootstrap smoke test를 통과했습니다. 공개 소스는 `https://github.com/K1mYongrok/scalableCounterfactual`에서 관리합니다. 원 라이선스 Stata 교차검증, 다중 운영체제 CI 및 대규모 Monte Carlo coverage 검증은 향후 방법론적 검증 과제로 남아 있습니다.

상세 내용은 `inst/doc/RELEASE_VALIDATION_1.0.0.md`를 참고하십시오.
