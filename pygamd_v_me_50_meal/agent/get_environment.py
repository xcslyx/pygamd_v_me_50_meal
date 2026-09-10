import os
import platform
import subprocess

from semantic_kernel.functions import kernel_function


class GetEnvironment:
    def __init__(self):
        pass

    @kernel_function(
        name="check_environment",
        description="获取当前系统的环境信息、当前工作目录以及指定目录下的文件列表。当需要确认文件是否存在、或者不知道当前处于哪个文件夹时，调用此工具。",
    )
    def check_environment(self, target_dir: str = ".") -> str:
        """
        参数:
        :param target_dir: 需要查看的目录路径。默认为 "." 表示当前目录。
        """
        try:
            # 1. 获取系统信息
            sys_info = f"操作系统: {platform.system()} {platform.release()}"
            
            # 2. 获取当前绝对路径
            current_pwd = os.getcwd()
            abs_target_dir = os.path.abspath(target_dir)
            
            # 3. 获取目标目录下的文件列表 (使用 ls -la 获取详细信息)
            # 相当于在终端执行 ls -la
            ls_result = subprocess.run(["ls", "-la", abs_target_dir], capture_output=True,  text=True, check=True)
            files_info = ls_result.stdout
            
            # 将信息组装成一段完整的报告返回给大模型
            report = (
                f"=== 系统状态 ===\n{sys_info}\n\n"
                f"=== 路径信息 ===\n"
                f"当前 Agent 运行路径: {current_pwd}\n"
                f"目标查看路径: {abs_target_dir}\n\n"
                f"=== 目录 [{abs_target_dir}] 下的文件列表 ===\n"
                f"{files_info}"
            )
            return report
            
        except Exception as e:
            return f"获取环境信息失败。错误: {str(e)}"