

def get_slack_info():
    import getpass
    from slackclient import SlackClient
    # import os
    try:
        from slack_creds import slack_token
    except Exception:
        return None

    # slack_token = os.environ["SLACK_API_TOKEN"]
    sc = SlackClient(slack_token)
    user = getpass.getuser()
    if user.startswith('konrad'): user = 'konradjk'
    users = [x['name'] for x in sc.api_call("users.list")['members']]
    default_channel = '#gnomad' if user not in users else '@' + user
    return sc, default_channel


def get_slack_channel_id(sc, channel):
    channel_id = [x for x in sc.api_call('channels.list')['channels'] if x['name'] == channel]
    return channel_id[0]['id'] if len(channel_id) and 'id' in channel_id[0] else None


def get_slack_user_id(sc, user):
    channel_id = [x for x in sc.api_call('users.list')['members'] if x['name'] == user]
    return channel_id[0]['id'] if len(channel_id) and 'id' in channel_id[0] else None


def send_message(channel=None, message="Your job is done!", icon_emoji=':woohoo:'):
    sc, default_channel = get_slack_info()

    sc.api_call(
        "chat.postMessage",
        channel=channel,
        text=message,
        icon_emoji=icon_emoji,
        parse='full'
    )


def send_snippet(channel=None, content='', filename='data.txt'):
    sc, default_channel = get_slack_info()

    get_channel = channel if channel is not None else default_channel
    if get_channel.startswith('@'):
        channel_id = get_slack_user_id(sc, get_channel.lstrip('@'))
    else:
        channel_id = get_slack_channel_id(sc, get_channel.lstrip('#'))

    try:
        sc.api_call("files.upload",
                    channels=channel_id,
                    content=content,
                    filename=filename)
    except Exception:
        print 'Slack connection fail. Was going to send:'
        print content


def try_slack(target, func, *args):
    import sys
    import os
    import traceback
    import time
    process = os.path.basename(sys.argv[0])
    try:
        func(*args)
        send_message(target, 'Success! {} finished!'.format(process))
    except Exception as e:
        if 'SparkException' in e:
            send_snippet(target, traceback.format_exc(), filename='error_{}_{}.txt'.format(process, time.strftime("%Y-%m-%d_%H:%M")))
        else:
            send_message(target, 'Job ({}) failed :white_frowning_face:\n```{}```'.format(process, traceback.format_exc()), ':white_frowning_face:')
        raise e