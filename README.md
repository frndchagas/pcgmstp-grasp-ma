# GRASP-MA para o PC-GMSTP

Implementação independente do GRASP com memória adaptativa descrito no capítulo 4 da tese de **Bruno José da Silva Barros (UFF, 2025)**. O objetivo deste repositório é facilitar a revisão do algoritmo pelo autor.

Todo o programa está em [`grasp.cpp`](grasp.cpp): leitura de instâncias GTP, construção, busca local, memória adaptativa, redução de arestas, linha de comando e testes pequenos. Usa apenas a biblioteca padrão do C++17. Não depende de CPLEX, mineradores, VPN ou infraestrutura de experimentos.

**Este não é o código original do Bruno.** É uma extração da nossa reimplementação, com as interpretações descritas abaixo. A correspondência com a tese ainda precisa ser validada pelo autor.

## Compilar e experimentar

No macOS com as ferramentas de linha de comando do Xcode, ou no Linux com GCC/Clang:

```sh
c++ -std=c++17 -O3 -Wall -Wextra -Wpedantic grasp.cpp -o grasp
./grasp --self-test
./grasp --example > exemplo.GTP
./grasp exemplo.GTP --seed 1 --iterations 12 --trace
```

O exemplo é sintético, tem dois grupos e melhor objetivo **−15**. Para uma instância própria:

```sh
./grasp caminho/instancia.GTP --seed 123 --iterations 300 --cpu-limit 600
```

Por padrão, o programa executa `max(300, 3 × número_de_grupos)` iterações, com redução do Teorema 11 e busca local acelerada. Esse perfil segue a seção 4.7.1 da tese; o Algoritmo 8 apresenta o perfil anterior `max(200, 2k)`. O número de iterações pode ser escolhido explicitamente.

Opções para conferir as variantes:

```sh
./grasp exemplo.GTP --seed 1 --iterations 12 --full-mst --trace
./grasp exemplo.GTP --seed 1 --iterations 12 --no-reduction --trace
./grasp --help
```

`--full-mst` recalcula a MST de cada vizinho, em vez de usar a floresta parcial. `--no-reduction` mantém todas as arestas entre grupos. A construção e a regra da memória permanecem iguais.

## Algoritmo e escolhas para revisão

Cada solução escolhe exatamente um representante de cada grupo. Seu objetivo é minimizar **custo da MST dos representantes − soma dos prêmios dos representantes**. As distâncias são euclidianas arredondadas ao inteiro mais próximo.

| Parte | Comportamento implementado |
|---|---|
| Construção | Roleta independente em cada grupo: probabilidade proporcional ao peso do vértice. Pesos iniciais iguais a 1. Não há RCL ou parâmetro alpha de construção. |
| Busca local | Embaralha os grupos, procura a melhor troca estritamente melhor dentro de cada grupo e aplica essa troca antes de visitar o próximo. Repete passagens até nenhuma melhora. |
| Memória: melhora global estrita | Reinicia todos os pesos em 1; os representantes da nova melhor solução recebem `tamanho_do_grupo/2`. |
| Memória: sem melhora global | Incrementa em 1 o peso de cada representante da candidata após a busca local, inclusive quando o objetivo empata. |
| Meios e epsilon | Pesos inteiros dobrados preservam os meios: inicial 2, incremento 2, reset do selecionado igual ao tamanho do grupo. `epsilon = 0`. |
| Aceleração, Teorema 10 | Calcula a floresta mínima dos demais representantes e acrescenta todas as arestas retidas incidentes ao candidato. Uma floresta residual desconexa é tratada por componentes. |
| Redução, Teorema 11 | Um grupo intermediário como testemunha; exclusões sequenciais por custo decrescente e, em empate, endpoints crescentes. Uma aresta ausente invalida o grupo como testemunha para aquele endpoint. |
| Aleatoriedade | `mt19937_64`, amostragem inteira por rejeição sem viés e Fisher–Yates. Vértices de cada grupo são examinados por ID crescente. |
| Empates | Apenas melhora estrita altera a solução. Entre trocas igualmente boas, fica a primeira examinada. A melhor global conserva a primeira candidata com aquele objetivo. |

Há detalhes da tese que exigiram interpretação:

- Os nomes alpha/beta aparecem trocados entre a prosa da seção 4.5.1 e o Algoritmo 8. Implementamos os comportamentos concretos de incremento e reset indicados na tabela.
- O Algoritmo 7 mostra uma passagem e omite detalhes da substituição do candidato e do reinício de suas arestas. A repetição até ótimo local segue o texto da seção 4.5.
- O PRNG, as sementes, os desempates e a ordem sequencial da redução são convenções desta implementação. Não temos o programa original para compará-los.
- A tese informa uma implementação em C. Este programa usa C++17. Coincidência de sementes ou tempos entre as duas implementações não é esperada.

Não há mineração, padrões, path relinking, elite ou exportação de frequências neste programa. Os únicos pesos que influenciam a construção são os da memória adaptativa. O limite de CPU e o trace são recursos de execução e inspeção desta implementação.

## Entrada, saída e parada

O leitor aceita o formato GTP das instâncias originais do dataset abaixo. Os campos são separados por espaços ou quebras de linha:

1. Quantidade de vértices e pares de coordenadas `(x,y)`.
2. Quantidade de grupos; para cada grupo, sua cardinalidade e os IDs dos vértices.
3. Um centro por grupo, pertencente ao respectivo grupo.
4. Sentinela `-999`, um prêmio por vértice e sentinela final `-999`.

Os IDs externos começam em **1**. Os grupos devem formar uma partição completa e disjunta. Coordenadas e prêmios devem ser inteiros; coordenadas têm módulo máximo de `10^9`, prêmios são não negativos e há no máximo 10.000 vértices. O programa armazena uma matriz de distâncias de tamanho quadrático; o limite de vértices não é uma recomendação de capacidade.

A saída é JSON. Com `--trace`, cada iteração concluída produz uma linha com a candidata após a busca local, seguida pela linha final `type: result`. Esta contém a melhor solução, seus representantes e arestas em IDs de vértices a partir de 1, objetivo, pesos finais dobrados, parâmetros e término. Não há promessa de ótimo global: o método é uma heurística.

O limite padrão de CPU é 600 segundos. Ele conta desde antes da leitura e é **cooperativo**: os checkpoints verificam leitura, redução e busca; operações de leitura, ordenação e escrita não são interrompidas no meio. Não é um limite físico imposto pelo sistema operacional. O trace também consome tempo. Para comparar trajetórias, use número fixo de iterações e um teto de CPU suficiente.

| Retorno | Significado |
|---:|---|
| 0 | Todas as iterações pedidas terminaram; ou comando informativo/self-test concluído. |
| 2 | Parada por CPU. `complete` é falso; a melhor candidata de uma iteração concluída é preservada, se existir. Se a parada ocorrer antes da busca, `best` é `null` e `termination` é `cpu_limit_before_search`. |
| 1 | Entrada, opções ou execução inválidas; mensagem em stderr. |

Uma iteração parcial não atualiza a memória nem substitui a melhor solução concluída. O programa não escreve resultados por conta própria: redirecionamentos de stdout são escolhidos por quem o executa.

## Verificações

`--self-test` usa enumeração de árvores por sequências de Prüfer como oráculo independente da MST por Kruskal. Confere sete instâncias sintéticas, 142 seleções de representantes nos grafos original/reduzido, 672 ótimos locais concluídos e 11 entradas inválidas. Inclui prêmios que tornam o objetivo negativo, um único grupo, arestas de custo zero, empates, coordenadas extremas e floresta parcial desconexa. Também confere a memória e a igualdade das trajetórias com MST completa e acelerada.

A extração foi comparada no Mac com o núcleo de origem: **24 combinações** de entrada, semente, redução e modo de MST, com **360 iterações** no total. Todas as candidatas, árvores, melhores soluções, pesos finais e quantidades de arestas removidas coincidiram. Isso verifica a extração do nosso núcleo; não demonstra equivalência com a implementação original do Bruno nem compara desempenho com a tese.

A automação deste repositório compila e executa os testes em macOS e Linux, incluindo AddressSanitizer e UndefinedBehaviorSanitizer. As verificações locais não usam resultados da campanha experimental.

## Referências

- Barros, Bruno José da Silva. **Variações do Problema da Árvore Geradora Mínima: Abordagens Algorítmicas e Estudo de Complexidade**. UFF, 2025. [PDF oficial](https://www.ic.uff.br/wp-content/tesesedissertacoes/frontend-tesesdissertacoes/download.php?id=1242.pdf&tipo=trabalho). Capítulo 4: Teorema 10 e Algoritmo 7, pp. 65–68; Teorema 11, pp. 68–69; GRASP-MA e Algoritmo 8, pp. 71–73; perfil acelerado, pp. 81–82. Numeração impressa da tese.
- [Instâncias PC-GMSTP e descrição do formato — Mendeley Data, versão 1](https://data.mendeley.com/datasets/8k4zs74njz/1). Os dados têm licença CC BY 4.0. Essa licença se refere aos dados; não é uma licença para o código original do autor. O repositório só inclui um exemplo sintético gerado por `--example`.
