from semantic_kernel.functions import kernel_function


class PygamdAnalysis:
    def __init__(self):
        pass
    
    @kernel_function(
        name="process_coordinates",
        description="处理坐标数据, 移除离子, 移除凝聚体的 PBC 等",
    )
    def process_coordinates(self, path: str, remove_ions_zhy: bool = False, remove_condensate_pbc: bool = False, lang: str = "zh") -> dict:
        """
        Args:
            path: 系统目录路径
            remove_ions_zhy: 是否移除离子, 默认 False
            remove_condensate_pbc: 是否移除凝聚体的 PBC, 默认 False
            lang: 语言，默认 "zh"
        """
        from pygamd_v_me_50_meal.data import Data
        data = Data(path, lang=lang)

        try:
            from pygamd_v_me_50_meal.pygamd_analysis.coordinates_processor_agent import CoordinatesProcessorAgent
            CoordinatesProcessorAgent(path, data, remove_ions_zhy, lang=lang).cal_xyz(remove_condensate_pbc)
            return {"messages": ["成功处理坐标数据"]}
        except Exception as e:
            return {"messages": [f"处理坐标数据失败出错: {e}"]}
    
    @kernel_function(
        name="cal_contact_map",
        description="计算分子间接触图 Contact Map",
    )
    def cal_contact_map(self, path: str, gpu_choice: str | int=0, r_cut: float=4.0,
                        cm_class_list: str="", balance_cut: str | None=None, domain: str='',
                        draw_limit: bool=False):
        """
        Args:
            path: 系统目录路径
            gpu_choice: GPU 选择, 默认 0
            r_cut: 截止距离, 默认 4
            cm_class_list: 分子组合列表, 指要计算接触图的分子组合，如 "1-1,1-3,2-2"。默认 ""。
            balance_cut: 选择平衡后的轨迹进行计算，如 "1000,2000"。默认 None。
            domain: 要计算的结构域, 如 "159-522,463-600"。默认 ""
            draw_limit: 是否绘制限制, 默认 False
        
        注意：在分子组合选择中，分子标号从 1 开始。
        """
        from pygamd_v_me_50_meal.data import Data
        data = Data(path)

        import multiprocessing
        multiprocessing.set_start_method('spawn', force=True)

        try:
            from pygamd_v_me_50_meal.pygamd_analysis.contact_map_calculator_agent import ContactMapCalculatorAgent
            ContactMapCalculatorAgent(path, data, gpu_choice=gpu_choice, r_cut=r_cut,
                                cm_class_list=cm_class_list, balance_cut=balance_cut, domain=domain,
                                draw_limit=draw_limit).calculate_contact_map_parallel()
            return {"messages": ["成功计算 Contact Map"]}
        except Exception as e:
            return {"messages": [f"计算 Contact Map 失败出错: {e}"]}

    @kernel_function(
        name="cal_rg",
        description="计算分子回转半径 Rg",
    )
    def cal_rg(self, path: str, cal_class_rg: str = "1", balance_cut: str | None = None, domain: str | None = None, calculate_mass: bool = True):
        """
        Args:
            path: 系统目录路径
            cal_class_rg: 要计算的分子, 指要计算回转半径的分子组合，如 "1,2"。默认 "1"。
            balance_cut: 选择平衡后的轨迹进行计算，如 "1000,2000"。默认 None, 即全部计算。
            domain: 要计算的结构域, 如 "159-522"。默认 None。注: 只能传入单个结构域。
            calculate_mass: 计算回转半径时是否考虑分子质量。默认 True。
        """
        try:
            from pygamd_v_me_50_meal.data import Data
            data = Data(path)

            from pygamd_v_me_50_meal.pygamd_analysis.rg_calculator_agent import RgCalculatorAgent
            RgCalculatorAgent(path, data, cal_class_rg=cal_class_rg, balance_cut=balance_cut, domain=domain, calculate_mass=calculate_mass).calculate()
            return {"messages": ["成功计算 Rg"]}
        except Exception as e:
            return {"messages": [f"计算 Rg 失败出错: {e}"]}

    @kernel_function(
        name="analysis_molecular_topology",
        description="分析分子拓扑结构，可统计键长、键角、二面角信息。",
    )
    def analysis_molecular_topology(self, path: str, lang: str = "zh", bond: bool = False, angle: bool = False, dihedral: bool = False):
        """ 
        Args:
            path: 系统目录路径
            lang: 语言，默认 "zh"
            bond: 是否统计键长, 默认 False
            angle: 是否统计键角, 默认 False
            dihedral: 是否统计二面角, 默认 False
        """
        from pygamd_v_me_50_meal.data import Data
        data = Data(path)

        try:
            if bond:
                from pygamd_v_me_50_meal.pygamd_analysis.molecule_topology_analysys.bond_analysys import BondAnalysys
                BondAnalysys(path, data).get_bond_length_parallel()

            if angle:
                from pygamd_v_me_50_meal.pygamd_analysis.molecule_topology_analysys.angle_analysys import AngleAnalysys
                AngleAnalysys(path, data, lang).get_angle_degree_parallel()

            if dihedral:
                from pygamd_v_me_50_meal.pygamd_analysis.molecule_topology_analysys.dihedral_analysys import DihedralAnalysys
                DihedralAnalysys(path, data, lang).get_dihedral_degree_parallel()
        except Exception as e:
            return {"messages": [f"分析分子拓扑结构失败出错: {e}"] }

    @kernel_function(
        name="cal_rmsd",
        description="计算分子均方根位移 RMSD",
    )
    def cal_rmsd(self, path: str, ref: str, cal_class_rmsd: str='1', domain: str | None=None, balance_cut: str | None=None):
        """
        计算 RMSD
        Args:
            path: 系统目录路径
            ref: 参考结构, 必须提供。 
            cal_class_rmsd: 要计算的分子, 指要计算 RMSD 的分子组合，如 "1,2"。默认 "1"。
            domain: 要计算的结构域, 如 "159-522"。默认 None。注: 只能传入单个结构域。
            balance_cut: 选择平衡后的轨迹进行计算，如 "1000,2000"。默认 None, 即全部计算。
        """
        try:
            from pygamd_v_me_50_meal.pygamd_analysis.rmsd_calculator import RMSDCalculator
            rmsd_calculator = RMSDCalculator(path, ref, cal_class_rmsd, domain, balance_cut)
            rmsd_calculator.calculate()
            return {"messages": [f"成功计算 RMSD: {rmsd_calculator.rmsd_results}"]}
        except Exception as e:
            return {"messages": [f"计算 RMSD 失败出错: {e}"]}
