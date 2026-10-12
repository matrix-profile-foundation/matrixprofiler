# Relatório de Auditoria: Gaps Funcionais, Assimetrias e Contratos de Retorno
**Package:** `matrixprofiler` (Computational Engine)  
**Data:** 10 de Outubro de 2026  
**Status da Migração `RcppThread`:** Concluída com 100% de sucesso (TBB e RcppParallel erradicados).

---

## 1. Sumário Executivo

O pacote `matrixprofiler` foi completamente saneado de dependências legadas de `RcppParallel` e Intel TBB, migrando 100% das suas rotinas multithread para `RcppThread`.
- **Compilação e Linkagem:** Limpa com `-latomic -lpthread`, sem dependências externas de bibliotecas TBB ou headers do `rcpp_parallel_fix.h`.
- **Testes Unitários:** 1.238 testes executados e aprovados (`devtools::test()`), com zero erros e zero warnings.
- **R CMD Check:** `0 errors | 0 warnings | 1 note` (nota inofensiva do compilador Ubuntu do sistema).
- **Fork Safety e Concorrência:** Com a eliminação do TBB e validação POSIX sob processos forkeados (`mclapply`, `future::multicore`), o fallback seguro para 1 worker evita crashes de runtime do R.
- **Detecção de CUDA:** Função utilitária exportada `has_cuda()` implementada para permitir que o usuário e pacotes clientes (`tsmp`) consultem em tempo de execução se o binário foi compilado com aceleração CUDA (`-DMATRIXPROFILER_ENABLE_CUDA`).

---

## 2. Matriz de Paridade Funcional entre Algoritmos

A tabela a seguir consolida o estado de suporte de cada funcionalidade transversal entre os algoritmos disponíveis no pacote:

| Algoritmo | Self-Join (AA) | AB-Join | Multithread CPU (`RcppThread`) | CUDA GPU | Suporte a NA/NaN/Inf | Anytime (`s_size < 1`) | Convenção de Índice Inválido | Perfis Direcionais (L/R) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **MPX** | ✅ Sim | ✅ Sim | ✅ Sim (`parallelFor`) | ✅ Sim (`mpx_gpu`) | ✅ Sim (`na_strategy`) | ✅ Sim (Dual-Path) | `NA_integer_` | ❌ Não |
| **STOMP** | ✅ Sim | ✅ Sim | ✅ Sim (`RcppThread`) | ❌ Não | ✅ Sim (`skip_location`) | ⚠️ *Exato (Ctrl+C parcial)* | `NA_integer_` | ✅ Sim (`left_right_profile`) |
| **STAMP** | ✅ Sim | ✅ Sim | ✅ Sim (`ThreadPool`) | ❌ Não | ⚠️ *Parcial* | ✅ Sim (Anytime) | `NA_integer_` | ❌ Não |
| **SCRIMP** | ✅ Sim | ✅ Sim | ✅ Sim (`RcppThread`) | ❌ Não | ⚠️ *Parcial* | ✅ Sim (Anytime) | `NA_integer_` | ❌ Não |
| **MASS** | N/A | ✅ Sim | ⚠️ *Apenas v3 normalized* | ❌ Não | ❌ Não | ❌ Não | N/A | N/A |
| **CONTRAST** | N/A | ✅ Sim | ✅ Sim | ❌ Não | ❌ Não | ❌ Não | N/A | N/A |
| **RFCP** | N/A | ✅ Sim | ✅ Sim | ✅ Sim (`rfcp_gpu`) | ✅ Sim | ✅ Sim (Top-K / Chunked) | `NA_integer_` | N/A |
| **WINDOWFUNC** | N/A | N/A | ✅ Sim (`muinvn`) | ❌ Não | ❌ Não | ❌ Não | N/A | N/A |

---

## 3. Detalhamento dos Gaps e Assimetrias Identificadas

### 3.1 Assimetrias de Concorrência e Multithreading

1. **`scrimp()` AB-Join Multithread e Pre-SCRIMP [RESOLVIDO]**:
   - Implementado `scrimpab_rcpp_parallel()` consumindo `mpxab_na_rcpp_parallel()`, eliminando o fallback serial e permitindo paralelismo em AB-joins.
   - `scrimp_rcpp_parallel()` atualizado para receber e aplicar `pre_scrimp` com `apply_pre_scrimp()`.
   - Adicionado fallback defensivo para fork safety sob processos POSIX (`mclapply`).

2. **`mass()` Paralelo Restrito a `type = "normalized"` e `version = "v3"`**:
   - Em `R/mass.R`, `n_workers > 1` só é repassado para o backend C++ se `type == "normalized"` e `version == "v3"`.
   - Se o usuário selecionar `type = "weighted"`, `"non_normalized"` ou `"absolute"`, ou `version = "v2"`, a execução é sempre serial em uma única thread, mesmo com `n_workers = 16`.

---

### 3.2 Padronização de Índices Inválidos (`NA_integer_`) [RESOLVIDO]

- **Status**: ✅ **100% Padronizado**.
- **Resolução**: `stamp.cpp` e `stomp.cpp` foram atualizados para inicializar buffers seriais, paralelos e direcionais (`left_profile_index`, `right_profile_index`) com `NA_INTEGER` (`NA_integer_` no R). Adicionalmente, corrigiu-se o offset de `left_profile_index` e `right_profile_index` em STOMP para 1-based (`i + 1`).
- **Impacto no R**: Elimina qualquer risco de indexação negativa acidental (`data[-1]`), garantindo total uniformidade de contratos entre `mpx`, `scrimp`, `stomp`, `stamp` e `rfcp`.

---

### 3.3 Relaxamento da Restrição de Tamanho de Janela (`window_size`) [RESOLVIDO]

- **Status**: ✅ **100% Implementado**.
- **Resolução**: A verificação rígida legada `window_size > ceiling(query_size / 2)` foi substituída pelas restrições matemáticas reais do algoritmo:
  - **AB-Join**: Permite $w \le \min(|A|, |B|) - 1$, possibilitando a busca de padrões usando quase a totalidade da query.
  - **Self-Join**: Permite $w < |data|$, validando dinamicamente se sobram subsequências fora da zona de exclusão efetiva ($(N - w) > \text{round}(w \times \text{ez})$).
- **Cobertura**: Atualizado em `mpx.R`, `stomp.R`, `stamp.R`, `scrimp.R` e `contrast.R`, com testes unitários em `test-algorithms.R`.

---

### 3.4 Tratamento de Dados Ausentes (`NA`, `NaN`, `Inf`)

1. **STOMP com Suporte a NA e Inf [RESOLVIDO]**:
   - `src/stomp.cpp` atualizado para rastrear `skip_data` e `skip_query` independentemente.
   - Entradas clonadas de forma defensiva para evitar mutação indesejada dos vetores originais.
   - `R/stomp.R` flexibilizado para `checkmate::qassert(data, "n+")` e `query, "n>=4"`.
   - Janelas com dados ausentes produzem distância `Inf` e índice `NA_integer_`, idêntico ao contrato do `mpx`.

2. **Funções de Janela Móvel (`windowfunc.R`)**:
   - `mov_mean`, `mov_var`, `mov_sum`, `mov_std`: Não suportam `NA`. Se um único valor `NA` ou `NaN` entra na janela, acumuladores recursivos (como Ogita ou soma móvel) ficam contaminados para sempre até o final do vetor.

---

### 3.5 Separação de Papéis: Motor Base (`matrixprofiler`) vs Mineração (`tsmp`)

Atualmente, certas funções do `matrixprofiler` realizam tarefas de mineração e detecção que pertencem conceitualmente ao `tsmp`:
- `contrast()`: Calcula e retorna metadados de `plato`, `plato_nn`, `plato_idx`, `plato_nn_idx`.
- `rfcp()`: Seleciona e materializa diretamente as coordenadas do platô ótimo.
- **Proposta de Contrato Mínimo de Retorno**:
  Para atuar como motor de ultra-performance, a camada base deve garantir saídas enxutas e sem custos de materialização de dados desnecessários:
  ```r
  list(
    matrix_profile      = <numeric vector>,
    profile_index       = <integer vector>,   # com NA_integer_ para invalidos
    left_matrix_profile = <numeric vector ou NULL>,
    right_matrix_profile= <numeric vector ou NULL>,
    partial             = <logical scalar>
  )
  ```
  Toda a lógica de identificação de *motifs*, *discords*, *platos* e regras de classificação pode ser centralizada no `tsmp`.

---

## 4. Próximos Passos Recomendados

1. **Unificação de Índices**: Substituir sistematicamente o `-1L` por `NA_integer_` nas saídas de `stamp` e `stomp`.
2. **Harmonização do `window_size`**: Flexibilizar a restrição de $50\%$ para permitir $w \le N - 1$ quando aplicável.
3. **Completar AB-Join Paralelo em SCRIMP**: Aproveitar o `mpxab_na_rcpp_parallel` já existente para habilitar execução paralela direta.
4. **Validação de Estresse de Memória (Sanitizers)**: Manter rotinas de checagem com ASAN/UBSan habilitados em builds de integração contínua.
5. **Aceleração C++ de Mineração Pós-MP (Opção B - Futuro)**:
   - Implementar motores em C++ no `matrixprofiler` para transferir os gargalos que hoje residem em loops puros do `tsmp`:
     - `mdl_subsequences_cpp()`: Discretização e seleção gulosa por MDL para `salient_subsequences()`.
     - `golden_section_cpp()` / `sdts_train_cpp()`: Otimização unidimensional e busca de padrões combinatórios para SDTS.
6. **Exposição e Contrato de GPU para o `tsmp`**:
   - Garantir estabilidade dos contratos de retorno e verificação de hardware (`has_cuda()`, `backend = "cuda"` em `mpx` e `rfcp`) para viabilizar o despacho transparente pelo `tsmp` com fallback automático para CPU.
