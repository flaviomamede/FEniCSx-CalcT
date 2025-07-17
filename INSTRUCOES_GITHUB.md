# Instruções para Atualizar o GitHub

Para atualizar o repositório GitHub com todas as implementações realizadas, siga os passos abaixo:

## 1. Execute o script de atualização

```bash
./atualizar_repositorio.sh
```

Este script vai:
- Adicionar todos os arquivos novos e modificados ao Git
- Criar um commit com uma mensagem descritiva
- Preparar tudo para o push

## 2. Faça o push para o GitHub

Após executar o script, envie as alterações para o GitHub com:

```bash
git push origin main
```

Você precisará inserir suas credenciais do GitHub quando solicitado.

## 3. Conteúdo da atualização

Esta atualização contém:

1. **Novo diretório CalcT-AbacusAI/** com a implementação completa e robusta da construção em camadas
   - Classe `ConstrucaoCamadasEstavel` para simulação estável
   - Testes completos de funcionalidade
   - Documentação detalhada

2. **Atualizações no diretório CalcT/**
   - Melhorias nos scripts R1 e R2
   - Documentação atualizada
   - Logs de simulação

## 4. Verificação

Após o push, acesse seu repositório no GitHub para verificar se todas as alterações foram aplicadas corretamente.