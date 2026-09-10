import os
import httpx
import getpass
import asyncio

from openai import AsyncOpenAI
from typing import Annotated

from semantic_kernel import Kernel
from semantic_kernel.connectors.ai.open_ai import OpenAIChatCompletion
from semantic_kernel.connectors.ai.function_choice_behavior import FunctionChoiceBehavior
from semantic_kernel.connectors.ai.chat_completion_client_base import ChatCompletionClientBase
from semantic_kernel.contents.chat_history import ChatHistory
from semantic_kernel.functions.kernel_arguments import KernelArguments

from semantic_kernel.connectors.ai.open_ai.prompt_execution_settings.open_ai_prompt_execution_settings import (
    OpenAIChatPromptExecutionSettings,
)

DEEPSEEK_API_KEY = os.getenv("DEEPSEEK_API_KEY")
if not DEEPSEEK_API_KEY:
    DEEPSEEK_API_KEY = getpass.getpass("Enter your DeepSeek API key: ")

custom_http_client = httpx.AsyncClient(timeout=120.0)
openai_client = AsyncOpenAI(
        api_key=DEEPSEEK_API_KEY,
        base_url="https://api.deepseek.com/v1",
        http_client=custom_http_client,
        max_retries=2  # 恢复原代码中的重试机制
    )

# Initialize the kernel
kernel = Kernel()

# Add OpenAI chat completion
kernel.add_service(OpenAIChatCompletion(
    ai_model_id="deepseek-v4-flash",
    async_client=openai_client
))

import logging

# Set the logging level for  semantic_kernel.kernel to DEBUG.
logging.basicConfig(
    format="[%(asctime)s - %(name)s:%(lineno)d - %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logging.getLogger("kernel").setLevel(logging.DEBUG)

chat_completion : OpenAIChatCompletion = kernel.get_service(type=ChatCompletionClientBase)

# Add the plugin to the kernel
from pygamd_v_me_50_meal.agent.pygamd_analysis import PygamdAnalysis
kernel.add_plugin(PygamdAnalysis(), plugin_name="PygamdAnalysis")

execution_settings = OpenAIChatPromptExecutionSettings()
execution_settings.function_choice_behavior = FunctionChoiceBehavior.Auto()

system_prompt = """
    你是一个名为 pygamd_v_me_50_meal 的分子动力学模拟专家助手。你的任务是理解用户的自然语言需求，并调用相关工具执行建模和数据分析。

    核心工作流规则（必须严格遵守）：
    1. 依赖检查：在执行任何分析工具（如计算 Contact Map、Rg、RMSD) 之前, 必须确认是否已经处理了坐标数据。如果不能确定，则首先进行提取坐标。
    2. PBC 处理：如果用户要求计算“质量数密度分布”，在调用提取坐标工具时，必须将 `remove_condensate_pbc` 参数设为 True。
    3. 目录校验：当用户提供体系名称时，必须验证其是否符合 "数字+分子名称-分子长度" 格式（如 40A-256+20B-512)。如果不符合，拒绝执行并指导用户修正。
    4. 遇到不确定的情况，必须要求用户确认，而不是自己判断。

    以清晰、专业的语言向用户报告执行结果。如果工具报错，请分析错误日志并提供修改建议。
    """ 
    
history = ChatHistory()
history.add_system_message(system_prompt)

# 1. Wrap the execution in an async function
async def main():
    # ... setup your kernel, custom HTTP client, and chat_completion ...
    
    while True:
        user_input = input("请输入您的需求: ")
        if user_input.lower() in ["quit", "exit", "q", "bye"]:
            print("Goodbye!")
            break
        
        history.add_user_message(user_input)
        
        try:
            # 2. 'await' is now safely inside an 'async def' function
            result = await chat_completion.get_chat_message_content(
                chat_history=history,
                settings=execution_settings,
                kernel=kernel
            )
            print(result.content)
            # history.add_assistant_message(...)
            
        except Exception as e:
            print(f"助手执行失败: {e}")

# 3. Use asyncio.run() to trigger the async function
if __name__ == '__main__':
    asyncio.run(main())
