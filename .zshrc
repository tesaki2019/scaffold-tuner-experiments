# alias python='python3'
alias pip='pip3'


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/tesak/anaconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/tesak/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/tesak/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/tesak/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
export PS1="%n@%m %1~ %# "
# <<< conda initialize <<<


export PATH="$PATH:/Applications/moe2024.0604/bin"

export PATH="$PATH:/Applications/ccg/moe/moe2024.0604/bin"
