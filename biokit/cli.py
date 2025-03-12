# biokit/cli.py
# coding='utf-8'
# Author:Tang Hongzhen
# Email: tanghongzhen34@gmail.com


import argparse


def main():
    parser = argparse.ArgumentParser(description="BioKit 命令行工具")
    subparsers = parser.add_subparsers(dest="command", help="子命令")

    # 示例子命令：序列分析
    parser_seq = subparsers.add_parser("seq", help="序列分析")
    parser_seq.add_argument("-i", "--input", required=True, help="输入文件")
    parser_seq.add_argument("-o", "--output", help="输出文件")

    args = parser.parse_args()

    if args.command == "seq":
        print(f"处理序列: {args.input} -> {args.output}")


if __name__ == "__main__":
    main()
